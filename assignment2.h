#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <vector>

using namespace std;

#define ISSQUARE(x)  (sqrt(x) - floor(sqrt(x)) == 0)

typedef struct
{
    int id;
    double x;
    double y;
} City;

typedef struct
{
    double cost;
    vector<int> path;
} PathCost;

vector<vector<vector<City>>> distributeCities(int numCitiesPerBlock, int numBlocksInRow, int numBlocksInCol, int gridDimX, int gridDimY);
void distributeBlocks(vector<vector<vector<City>>> blockedCities);

double fRand(double fMin, double fMax)
{
    // https://stackoverflow.com/a/2704552
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void printMatrix(double **matrix, int r, int c)
{
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void printBlocked(vector<vector<vector<City>>> blocks)
{
    for (int i = 0; i < (int)blocks.size(); i++)
    {
        // printf("We're in block %i\n", i);
        printf("Block %i {\n", i);
        for (int j = 0; j < (int)blocks[i].size(); j++)
        {
            // printf("In row %i of block %i we have\n", j, i);
            printf("\t[");
            for (int k = 0; k < (int)blocks[i][j].size(); k++)
            {
                printf("%i:(%.2f, %.2f) ", blocks[i][j][k].id, blocks[i][j][k].x, blocks[i][j][k].y);
            }
            printf("]\n");
        }
        printf("}\n\n");
    }
}

void printMatrixArray(vector<City> matrix, int rowWidth, int numElements)
{
    int counter = 0;
    for (int i = 0; i < (numElements / (float)rowWidth); i++)
    {
        printf("[ ");
        for (int j = 0; j < min((numElements - (i * rowWidth)), rowWidth); j++)
        {
            City city = matrix[counter];
            printf("(%f, %f) ", city.x, city.y);
            counter++;
        }
        printf("]\n");
    }
}

typedef struct
{
    int blockId;
    vector<City> path;
    double cost;
} BlockSolution;

typedef struct
{
    int threadId;
    vector<City> cities;
} TSPArgs;

struct sortByX
{
    inline bool operator()(const City city1, const City &city2)
    {
        return (city1.x < city2.x);
    }
};

struct sortByY
{
    inline bool operator()(const City city1, const City &city2)
    {
        return (city1.y < city2.y);
    }
};

double distance(City c1, City c2)
{
    return sqrt(pow(c1.x - c2.x, 2) + pow(c1.y - c2.y, 2));
}

vector<City> readCities(char *filePath)
{
    vector<City> cities;
    ifstream file(filePath);
    double x, y;
    int id = 0;
    // Reading in the cities tried to methodize this and failed
    while (file >> x >> y)
    {
        City city;
        city.id = id;
        city.x = x;
        city.y = y;
        cities.push_back(city);
        id++;
    }
    printf("We just read in %lu cities\n", cities.size());
    return cities;
}

void genKey(vector<int> set, int z, long long &key)
{
    key = 0;
    key |= z;
    for (int j : set)
    {
        key |= (1 << (j + 8));
    }
}

vector<vector<int>> generateSubsets(int size, int n)
{
    int count = 0;
    vector<vector<int>> container;
    vector<int> row;
    vector<bool> v((unsigned long)n);
    fill(v.begin(), v.begin() + size, true);

    do
    {
        for (int i = 0; i < n; ++i)
        {
            if (v[i])
            {
                count++;
                row.push_back(i + 1);
                if (count == size)
                {
                    container.push_back(row);
                    row.clear();
                    count = 0;
                }
            }
        }
    } while (prev_permutation(v.begin(), v.end()));
    return container;
}

double **computeDistanceMatrix(vector<City> cities)
{
    double **distances = (double **)malloc((int)cities.size() * sizeof(double *));
    for (int i = 0; i < (int)cities.size(); i++)
        distances[i] = (double *)malloc((int)cities.size() * sizeof(double));

    for (int i = 0; i < (int)cities.size(); i++)
    {
        for (int j = 0; j < (int)cities.size(); j++)
        {
            City city1 = cities[i];
            City city2 = cities[j];
            distances[i][j] = sqrt(pow(city1.x - city2.x, 2) + pow(city1.y - city2.y, 2));
        }
    }
    return distances;
}

void printPath(vector<int> path)
{

    printf("path is: ");
    for (int i = 0; i < (int)path.size() - 1; i++)
    {
        printf("%i -> ", path[i]);
    }
    printf("%i", path[(int)path.size() - 1]);

    printf("\n");
}
