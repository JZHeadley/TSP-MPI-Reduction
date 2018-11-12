#include <string.h>
#include <time.h>
#include <limits.h>
#include <utility>
#include <cmath>

#include "assignment2.h"

#define BLOCK_NUM_TAG 6
#define NUM_BLOCKS_TAG 66
#define BLOCK_CITIES_TAG 666
#define NUM_CITIES_TAG 6666
#define CITIES_RESULT_TAG 66666
#define COST_RESULT_TAG 7
#define NUM_BLOCKS_RECV_TAG 77
#define PATH_COST_TAG 777

int numProcs;
int procNum;
int coords[2];
MPI_Datatype mpi_city_type;
MPI_Datatype mpi_block_solution_type;

void initMPI()
{
    const int numItems = 3;
    int blockLengths[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[3];
    offsets[0] = offsetof(City, id);
    offsets[1] = offsetof(City, x);
    offsets[2] = offsetof(City, y);
    MPI_Type_create_struct(numItems, blockLengths, offsets, types, &mpi_city_type);
    MPI_Type_commit(&mpi_city_type);

    MPI_Comm_rank(MPI_COMM_WORLD, &procNum);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
}

double pathDistance(vector<City> path)
{
    City previous = path[0];
    double pathCost = 0;
    for (int i = 1; i < path.size(); i++)
    {
        pathCost += distance(previous, path[i]);
        previous = path[i];
    }
    return pathCost;
}

BlockSolution MPI_ManualReduce(BlockSolution solution, MPI_Comm comm)
{
    // I came up with the same idea as this guy and am just using his code as a base and modifying it to my need
    // I just don't feel like trying to think about all the reduction math from scratch https://gist.github.com/rmcgibbo/7178576
    int tag = 0;
    int size;
    int rank;
    int recvNumCities;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    const int lastpower = 1 << (int)log2(size);
    MPI_Status status;
    City *recvCities;
    City *sendCities;
    double recvPathCost;
    vector<City> path;

    // each of the ranks greater than the last power of 2 less than size
    // need to downshift their data, since the binary tree reduction below
    // only works when N is a power of two.
    for (int i = lastpower; i < size; i++)
        if (rank == i)
        {
            // printf("sending blocks in %i to %i\n", rank, i - lastpower);
            int numCities = solution.path.size();
            MPI_Send(&numCities, 1, MPI_INT, i - lastpower, NUM_CITIES_TAG, comm);
            sendCities = (City *)calloc(sizeof(City), numCities);
            copy(solution.path.begin(), solution.path.end(), sendCities);
            MPI_Send(sendCities, numCities, mpi_city_type, i - lastpower, tag, comm);
            MPI_Send(&solution.cost, 1, MPI_DOUBLE, i - lastpower, PATH_COST_TAG, comm);
        }
    for (int i = 0; i < size - lastpower; i++)
    {
        if (rank == i)
        {
            MPI_Recv(&recvNumCities, 1, MPI_INT, i + lastpower, NUM_CITIES_TAG, comm, &status);
            printf("process %i is about to receive %i cities from process %i\n", rank, recvNumCities, i + lastpower);
            recvCities = (City *)calloc(sizeof(City), recvNumCities);
            MPI_Recv(recvCities, recvNumCities, mpi_city_type, i + lastpower, tag, comm, &status);
            MPI_Recv(&recvPathCost, 1, MPI_DOUBLE, i + lastpower, PATH_COST_TAG, comm, &status);
            BlockSolution recvBlock;
            for (int x = 0; x < recvNumCities; x++)
                path.push_back(recvCities[x]);
            recvBlock.path = path;
            recvBlock.cost = recvPathCost;
            recvBlock.blockId = rank;
            solution = mergeBlocks(solution, recvBlock);
        }
    }
    // here we do a proper reduction using powers of 2
    for (int d = 0; d < (int)log2(lastpower); d++)
        for (int k = 0; k < lastpower; k += 1 << (d + 1))
        {
            const int receiver = k;
            const int sender = k + (1 << d);
            if (rank == receiver)
            {
                MPI_Recv(&recvNumCities, 1, MPI_INT, sender, NUM_CITIES_TAG, comm, &status);
                recvCities = (City *)calloc(sizeof(City), recvNumCities);
                MPI_Recv(recvCities, recvNumCities, mpi_city_type, sender, tag, comm, &status);
                MPI_Recv(&recvPathCost, 1, MPI_DOUBLE, sender, PATH_COST_TAG, comm, &status);

                BlockSolution recvBlock;
                for (int x = 0; x < recvNumCities; x++)
                    path.push_back(recvCities[x]);
                recvBlock.path = path;
                recvBlock.cost = recvPathCost;
                recvBlock.blockId = rank;
                solution = mergeBlocks(solution, recvBlock);
            }
            else if (rank == sender)
            {
                int numCities = solution.path.size();
                MPI_Send(&numCities, 1, MPI_INT, receiver, NUM_CITIES_TAG, comm);
                sendCities = (City *)calloc(sizeof(City), numCities);
                copy(solution.path.begin(), solution.path.end(), sendCities);
                MPI_Send(sendCities, numCities, mpi_city_type, receiver, tag, comm);
                // MPI_Send(&solution.path[0], numCities, mpi_city_type, receiver, tag, comm);
                MPI_Send(&solution.cost, 1, MPI_DOUBLE, receiver, PATH_COST_TAG, comm);
            }
        }
    return solution;
}

vector<int> getBlocksPerDim(int numBlocks)
{
    int numBlocksInRow;
    int numBlocksInCol;
    if (ISSQUARE(numBlocks))
    {
        numBlocksInRow = sqrt(numBlocks);
        numBlocksInCol = sqrt(numBlocks);
    }
    else
    {
        int divisor = 2;
        while (numBlocks % divisor != 0)
        {
            // will either take no steps if we have an even number, 1 step if we have an odd, or n-2 steps if we have a prime number
            divisor++;
        }
        numBlocksInRow = divisor;
        numBlocksInCol = numBlocks / numBlocksInRow;
    }
    return vector<int>({numBlocksInRow, numBlocksInCol});
}

vector<vector<City>> distributeBlocks(vector<vector<City>> blocks, int numBlocks, int numCitiesPerBlock, MPI_Comm comm)
{
    int blocksLeft = numBlocks;
    int blockToSend = 0;
    MPI_Request request;
    int *blocksToSend = (int *)calloc(sizeof(int), numProcs);
    vector<vector<City>> leftovers;

    while (blocksLeft)
    {
        blocksToSend[blocksLeft % numProcs]++;
        blocksLeft--;
    }

    for (int i = 0; i < numProcs; i++) // for each worker process
    {
        // printf("sending %i blocks to process %i\n\n", blocksToSend[i], i);
        if (i > 0)
            MPI_Send(&blocksToSend[0] + i, 1, MPI_INT, i, NUM_BLOCKS_RECV_TAG, comm);
        for (int j = 0; j < blocksToSend[i]; j++)
        {
            if (i == 0)
            {
                leftovers.push_back(blocks[blockToSend]);
            }
            else
            {
                City *block = (City *)calloc(numCitiesPerBlock, sizeof(City));
                copy(blocks[blockToSend].begin(), blocks[blockToSend].end(), block);
                MPI_Send(block, numCitiesPerBlock, mpi_city_type, i, BLOCK_CITIES_TAG, comm);
            }
            blockToSend++;
        }
    }
    free(blocksToSend);
    return leftovers;
}

inline double swapPairCost(pair<City, City> left, pair<City, City> right)
{
    return (distance(left.first, right.second) + distance(left.second, right.first) - distance(left.first, left.second) - distance(right.first, right.second));
}

BlockSolution mergeBlocks(BlockSolution solution1, BlockSolution solution2)
{
    double bestSwapCost = INT_MAX;
    vector<City> cities1 = solution1.path;
    vector<City> cities2 = solution2.path;
    pair<City, City> cityPair1;
    pair<City, City> cityPair2;
    pair<pair<City, City>, pair<City, City>> bestSwapEdges;
    double swapCost;

    for (int i = 0; i < cities1.size(); i++)
    {
        cityPair1 = make_pair(cities1[0], cities1[1]);
        for (int j = 0; j < cities2.size(); j++)
        {
            cityPair2 = make_pair(cities2[0], cities2[1]);
            swapCost = swapPairCost(cityPair1, cityPair2);
            if (swapCost < bestSwapCost)
            {
                bestSwapCost = swapCost;
                bestSwapEdges = make_pair(cityPair1, cityPair2);
            }
            rotate(cities2.begin(), cities2.begin() + 1, cities2.end());
        }
        rotate(cities1.begin(), cities1.begin() + 1, cities1.end());
    }

    cities1 = solution1.path;
    cities2 = solution2.path;
    cities2.pop_back(); // remove the final node so we're not circular anymore
    // printPath(cities1);
    // printPath(cities2);
    // printf("best swap is from cities1 %i to cities2 %i and cities2 %i to cities1 %i\n",
        //    bestSwapEdges.first.first.id, bestSwapEdges.second.first.id, bestSwapEdges.second.second.id, bestSwapEdges.first.second.id);
    while (cities2[0].id != bestSwapEdges.second.first.id)
    {
        rotate(cities2.begin(), cities2.begin() + 1, cities2.end());
    }
    // do a final rotation so that the path to add is now essentially in order but backwards
    rotate(cities2.begin(), cities2.begin() + 1, cities2.end());

    // printPath(cities2);
    vector<City> path;
    bool flag = true;
    for (int i = 0; i < cities1.size(); i++)
    {
        path.push_back(cities1[i]);
        // we should be at the start of the new edge in left block
        // we also should already have block 2 rotated so the right half of this edge is at the head
        if ((cities1[i].id == bestSwapEdges.first.first.id || cities1[i].id == bestSwapEdges.first.second.id) && flag)
        {
            flag = false;
            for (int j = cities2.size() - 1; j >= 0; j--)
            {
                path.push_back(cities2[j]);
            }
        }
    }

    BlockSolution merged;
    merged.blockId = procNum;
    merged.cost = solution1.cost + solution2.cost + bestSwapCost; // - distance(bestSwapEdges.first.first, bestSwapEdges.first.second) - distance(bestSwapEdges.second.first, bestSwapEdges.second.second);
    merged.path = path;
    // printf("path size is %i\n", path.size());
    // printPath(path);
    // printf("\n");
    return merged;
}
int main(int argc, char **argv)
{
    time_t t;
    srand(0);

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    MPI_Init(&argc, &argv);
    initMPI();
    if (argc != 5)
    {
        printf("Usage:  ./tsp numCitiesPerBlock numBlocks gridDimX gridDimY\n");
        exit(1);
    }
    int numCitiesPerBlock = atoi(argv[1]);
    int numBlocks = atoi(argv[2]);
    int gridDimX = atoi(argv[3]);
    int gridDimY = atoi(argv[4]);
    if (numCitiesPerBlock > 16)
    {
        if (procNum == 0)
            printf("Come on... We don't want to wait forever so lets just have you retry that with less than 16 cities per block...\n");
        MPI_Finalize();
        exit(1337);
    }

    MPI_Comm grid_comm;
    vector<int> procDims = getBlocksPerDim(numProcs);
    int wrap[2] = {1, 1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, &procDims[0], wrap, 0, &grid_comm);
    MPI_Comm_rank(grid_comm, &procNum);
    MPI_Status status;

    MPI_Cart_coords(grid_comm, procNum, 2, coords);
    // printf("proc %i is at ( %i, %i )\n", procNum, coords[0], coords[1]);
    if (procNum == 0)
        printf("We have %i cities for each of our %i blocks\n", numCitiesPerBlock, numBlocks);
    vector<BlockSolution> blockSolutions;
    MPI_Request request;
    if (procNum == 0)
    {
        vector<int> blockedDims = getBlocksPerDim(numBlocks);
        // generate our cities in their respective blocks
        vector<vector<City>> blockedCities = distributeCities(numCitiesPerBlock, blockedDims[0], blockedDims[1], gridDimX, gridDimY);
        // Distribute our blocks to each processor saving some of the blocks for our main block
        vector<vector<City>> blocksForHead = distributeBlocks(blockedCities, numBlocks, numCitiesPerBlock, grid_comm);
        // deal with the leftover block
        for (int i = 0; i < blocksForHead.size(); i++)
        {
            blockSolutions.push_back(tsp(blocksForHead[i]));
        }
    }
    else
    {
        // Kill off uneccesary blocks
        if (procNum > numBlocks)
        {
            MPI_Finalize();
            return 0;
        }
        // receive our blocks for head process here and run tsp on them
        int numBlocksToRecv;
        MPI_Recv(&numBlocksToRecv, 1, MPI_INT, 0, NUM_BLOCKS_RECV_TAG, grid_comm, &status);
        while (numBlocksToRecv)
        {
            City *block = (City *)calloc(numCitiesPerBlock, sizeof(City));
            MPI_Recv(block, numCitiesPerBlock, mpi_city_type, 0, BLOCK_CITIES_TAG, grid_comm, &status);
            vector<City> cities;
            for (int x = 0; x < numCitiesPerBlock; x++)
            {
                cities.push_back(block[x]);
            }
            blockSolutions.push_back(tsp(cities));
            numBlocksToRecv--;
        }
    }
    // perform a local reduction within the process since we might have multiple blocks per process
    while (blockSolutions.size() > 1)
    {
        blockSolutions[0] = mergeBlocks(blockSolutions[0], blockSolutions[1]);
        blockSolutions.erase(blockSolutions.begin() + 1);
    }
    // perform the full reduction here
    MPI_Barrier(grid_comm);
    BlockSolution finalSolution = MPI_ManualReduce(blockSolutions[0], grid_comm);
    MPI_Barrier(grid_comm);

    if (procNum == 0)
    {
        clock_gettime(CLOCK_MONOTONIC_RAW, &end);
        uint64_t diff = (1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec) / 1e6;

        printf("TSP ran in %llu ms for %lu cities and the trip cost %f\n", (long long unsigned int)diff, numBlocks * numCitiesPerBlock, finalSolution.cost);
    }

    MPI_Finalize();
    return 0;
}

/**
 * Generates and distributes cities to the block structure based on the input parameters
 */
vector<vector<City>> distributeCities(int numCitiesPerBlock, int numBlocksInRow, int numBlocksInCol, int gridDimX, int gridDimY)
{
    vector<vector<City>> blockedCities;
    int numBlocks = numBlocksInRow * numBlocksInCol;
    printf("%i blocks in X %i in Y\n", numBlocksInRow, numBlocksInCol);
    float xSpacePerBlock = gridDimX / (float)numBlocksInRow;
    float ySpacePerBlock = gridDimY / (float)numBlocksInCol;
    int id = 0;
    City city;
    vector<City> cities;

    for (int i = 0; i < numBlocks; i++)
    {
        cities.clear();
        for (int j = 0; j < numCitiesPerBlock; j++)
        {

            city.id = id;
            int row = ((i - (i % numBlocksInRow)) / numBlocksInRow);
            int col;
            col = (numBlocksInCol - (i % numBlocksInCol)) - 1;
            city.x = fRand(row * xSpacePerBlock, (row + 1) * xSpacePerBlock);
            city.y = fRand(col * ySpacePerBlock, (col + 1) * ySpacePerBlock);

            cities.push_back(city);
            id++;
        }
        blockedCities.push_back(cities);
    }
    return blockedCities;
}

BlockSolution tsp(vector<City> cities)
{

    double **distances = computeDistanceMatrix(cities);
    map<long long int, PathCost> solutionsMap;
    vector<int> minPath{};
    double minCost = INT_MAX;

    int start = 0;
    int numCities = cities.size();
    long long key = 0x00000;
    vector<int> cityNums;
    // convert cities back to integer array
    for (int i = 1; i < numCities; i++)
    {
        cityNums.push_back(i);
    }

    // initalize first 2 levels of the lookup table
    for (int i = 1; i < numCities; i++)
    {
        for (int j = 1; j < numCities; j++)
        {
            if (i == j)
                continue;
            vector<int> iSet{i};
            genKey(iSet, j, key);
            PathCost pathCost;
            vector<int> path{0, i};
            pathCost.path = path;
            pathCost.cost = distances[i][j] + distances[0][i];
            solutionsMap.insert(pair<long long, PathCost>(key, pathCost));
        }
    }
    //we're good so far...
    double currentCost = 0;

    for (int i = 2; i < numCities; i++)
    { // iterate through all cardinalities of subsets
        // printf("working on subsets of size %i\n", i);
        vector<vector<int>> subsets = generateSubsets(i, cityNums.size());
        for (vector<int> set : subsets)
        {
            for (int k : set)
            {
                vector<int> kSet{k};
                vector<int> diff;
                set_difference(set.begin(), set.end(), kSet.begin(), kSet.end(), inserter(diff, diff.begin()));
                double minCost = INT_MAX;
                vector<int> minPath;
                int bestM;
                // we initialized 2 levels earlier so this for loop will always be able to run.
                for (int m : diff)
                {
                    vector<int> mSet{m}; // need to generate the key for k-1
                    vector<int> noMoreM; // get rid of m because thats where we're going
                    set_difference(diff.begin(), diff.end(), mSet.begin(), mSet.end(), inserter(noMoreM, noMoreM.begin()));

                    genKey(noMoreM, m, key);
                    currentCost = solutionsMap[key].cost + distances[m][k];
                    if (currentCost < minCost)
                    {
                        minCost = currentCost;
                        minPath = solutionsMap[key].path;
                        bestM = m;
                    }
                }
                genKey(diff, k, key);

                PathCost pathCost;
                pathCost.cost = minCost;
                minPath.push_back(bestM);
                pathCost.path = minPath;
                solutionsMap.insert(pair<long long, PathCost>(key, pathCost));
            }
        }
    }

    int bestM;
    for (int m : cityNums)
    {
        vector<int> mSet{m}; // need to generate the key for k-1
        vector<int> noMoreM; // get rid of m because thats where we're going
        set_difference(cityNums.begin(), cityNums.end(), mSet.begin(), mSet.end(), inserter(noMoreM, noMoreM.begin()));

        genKey(noMoreM, m, key);
        currentCost = solutionsMap[key].cost + distances[m][0];
        if (currentCost < minCost)
        {
            minCost = currentCost;
            vector<int> path = solutionsMap[key].path;
            minPath = path;
            bestM = m;
        }
    }
    free(distances);
    minPath.push_back(bestM);
    minPath.push_back(0);
    BlockSolution blockSolution;
    vector<City> truePath = convPathToCityPath(cities, minPath);
    blockSolution.path = truePath;
    blockSolution.cost = minCost;
    blockSolution.blockId = procNum;
    return blockSolution;
}