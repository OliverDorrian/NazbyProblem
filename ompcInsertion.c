#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "coordReader.h"

double calculateDistance(double x1, double y1, double x2, double y2) {
    return sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}

void createDistanceMatrix(double **passedCoords, int numCoordinates, double **distanceMatrix) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < numCoordinates; ++i) {
        for (int j = i + 1; j < numCoordinates; ++j) {
            double distance = calculateDistance(passedCoords[i][0], passedCoords[i][1], passedCoords[j][0],
                                                passedCoords[j][1]);
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
        }
    }

    for (int i = 0; i < numCoordinates; ++i) {
        distanceMatrix[i][i] = 0.000000;
    }
}

void cheapestInsertion(double **distanceMatrix, int numCoordinates, int *tour) {
    int unvisited[numCoordinates];
    for (int i = 0; i < numCoordinates; i++) {
        unvisited[i] = i;
    }

    int tourSize = 1;
    int startCity = 0;
    tour[0] = startCity;
    unvisited[startCity] = -1;  // Mark the starting city as visited

    for (int i = 1; i < numCoordinates; i++) {
        int bestCity = -1;
        int bestInsertionIndex = -1;
        double minCost = 1.0 / 0.0;

        #pragma omp parallel for schedule(static)
        for (int vk = 0; vk < numCoordinates; vk++) {
            if (unvisited[vk] != -1) {  // Check for unvisited cities
                double local_minCost = 1.0 / 0.0; // Equivalent to std::numeric_limits<double>::max() in C++
                int local_bestCity = -1;
                int local_bestInsertionIndex = -1;

                for (int vn = 0; vn < tourSize; vn++) {
                    int vn_1 = tour[vn];
                    int vn_2 = tour[(vn + 1) % tourSize];  // Circular tour

                    double insertionCost = distanceMatrix[vn_1][vk] + distanceMatrix[vn_2][vk] - distanceMatrix[vn_1][vn_2];

                    if (insertionCost < local_minCost) {
                        local_minCost = insertionCost;
                        local_bestCity = vk;
                        local_bestInsertionIndex = vn;
                    }
                }

                #pragma omp critical
                {
                    if (local_minCost < minCost) {
                       minCost = local_minCost;
                       bestCity = local_bestCity;
                       bestInsertionIndex = local_bestInsertionIndex;
                    }
                }

            }
        }

        tourSize++;  // Increment tour size
        int nextTourIndex = (bestInsertionIndex + 1) % tourSize;  // Circular tour

        // Move cities in the tour to make space for the bestCity
        for (int vn = tourSize - 1; vn > nextTourIndex; vn--) {
            tour[vn] = tour[vn - 1];
        }

        tour[nextTourIndex] = bestCity;
        unvisited[bestCity] = -1;  // Mark the chosen city as visited
    }
}

int main(int argc, char *argv[]) {

    double start, end, time_taken;

    int numThreads = 12;
    omp_set_num_threads(numThreads);

    argv[1] = "4096_coords.coord";
    char const *fileName = argv[1];

    argv[2] = "compOut.dat";
    //char const *fileName = "4096_coords.coord";
    //char const *fileName = "9_coords.coord";

    int numCoordinates = readNumOfCoords(fileName);
    double **coords = readCoords(fileName, numCoordinates);

    start = omp_get_wtime();

    // Allocate memory for the distance matrix
    double **distanceMatrix = (double **)malloc(numCoordinates * sizeof(double *));
    for (int i = 0; i < numCoordinates; i++) {
        distanceMatrix[i] = (double *)malloc(numCoordinates * sizeof(double));
    }

    // Generate distance Matrix
    createDistanceMatrix(coords, numCoordinates, distanceMatrix);

    int tour[numCoordinates];
    cheapestInsertion(distanceMatrix, numCoordinates, tour);

    end = omp_get_wtime();
    time_taken = ((double) (end - start));
    printf("Execution time: %.6f seconds\n", time_taken);

    writeTourToFile(tour, numCoordinates, argv[2]);

    // Free Memory
    for (int i = 0; i < numCoordinates; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);
    }

    free(coords);
    free(distanceMatrix);

    return 0;
}
