#include <string.h>
#include <mpi.h>
#include <time.h>

#include "assignment2.h"

int numProcs;
int procNum;
int coords[2];

void initMPI()
{
    const int numItems = 3;
    int blockLengths[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Datatype mpi_city_type;
    MPI_Aint offsets[3];
    offsets[0] = offsetof(City, id);
    offsets[1] = offsetof(City, x);
    offsets[2] = offsetof(City, y);
    MPI_Type_create_struct(numItems, blockLengths, offsets, types, &mpi_city_type);
    MPI_Type_commit(&mpi_city_type);

    MPI_Comm_rank(MPI_COMM_WORLD, &procNum);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
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

int main(int argc, char **argv)
{
    time_t t;
    srand(t);

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
    else if (numBlocks < numProcs)
    {
        if (procNum == 0)
            printf("Either run it with fewer processors or more blocks...I don't feel like writing the logic to kill excess processes\n");
        MPI_Finalize();
        exit(1477);
    }
    MPI_Comm grid_comm;
    vector<int> procDims = getBlocksPerDim(numProcs);
    int wrap[2] = {1, 1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, &procDims[0], wrap, 0, &grid_comm);
    MPI_Comm_rank(grid_comm, &procNum);

    MPI_Cart_coords(grid_comm, procNum, 2, coords);
    // printf("proc %i is at ( %i, %i )\n", procNum, coords[0], coords[1]);
    if (procNum == 0)
        printf("We have %i cities for each of our %i blocks\n", numCitiesPerBlock, numBlocks);

    if (procNum == 0)
    {
        vector<int> blockedDims = getBlocksPerDim(numBlocks);
        // generate our cities in their respective blocks
        vector<vector<vector<City>>> blockedCities = distributeCities(numCitiesPerBlock, blockedDims[0], blockedDims[1], gridDimX, gridDimY);
        // Distribute our blocks to each processor saving one of the blocks for our main block
        // distributeBlocks(blockedCities);

        // deal with the leftover block

        // figure out how to do a reduction across the rows using the TSP-merge technique
        // then do it across the columns

        // should have the result collected here in process 0 due to how reductions work
    }
    else
    {
        // receive our blocks for all other processes here and run tsp on them

        // might have to write the code to do this TSP-merge reduction ourselves...
    }
    if (procNum == 0)
    {
        clock_gettime(CLOCK_MONOTONIC_RAW, &end);
        uint64_t diff = (1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec) / 1e6;

        printf("TSP ran in %llu ms for %lu cities\n", (long long unsigned int)diff, numBlocks * numCitiesPerBlock);
    }

    MPI_Finalize();
    return 0;
}

/**
 * Generates and distributes cities to the block structure based on the input parameters
 */
vector<vector<vector<City>>> distributeCities(int numCitiesPerBlock, int numBlocksInRow, int numBlocksInCol, int gridDimX, int gridDimY)
{
    vector<vector<vector<City>>> blockedCities;

    printf("%i blocks in X %i in Y\n", numBlocksInRow, numBlocksInCol);
    float xSpacePerBlock = gridDimX / (float)numBlocksInRow;
    float ySpacePerBlock = gridDimY / (float)numBlocksInCol;
    int id = 0;
    for (int i = 0; i < numBlocksInRow; i++)
    {
        blockedCities.push_back(vector<vector<City>>());
        for (int j = 0; j < numBlocksInCol; j++)
        {
            blockedCities[i].push_back(vector<City>());
            for (int k = 0; k < numCitiesPerBlock; k++)
            {
                City city;
                city.id = id;
                city.x = fRand(i * xSpacePerBlock, (i + 1) * xSpacePerBlock);
                city.y = fRand(j * ySpacePerBlock, (j + 1) * ySpacePerBlock);
                // printf("city at (%f, %f) which should be between (%f, %f) and (%f, %f)\n",city.x, city.y, i * xSpacePerBlock, j * ySpacePerBlock, (i + 1) * xSpacePerBlock, (j + 1) * ySpacePerBlock);
                blockedCities[i][j].push_back(city);
                id++;
            }
        }
    }
    return blockedCities;
}
