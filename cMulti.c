#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>

#include "coordReader.c"

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
        double minCost = 1.0 / 0.0; // Equivalent to std::numeric_limits<double>::max() in C++

        #pragma omp parallel for schedule(static)
        for (int vk = 0; vk < numCoordinates; vk++) {
            if (unvisited[vk] != -1) {  // Check for unvisited cities
                double local_minCost = 1.0 / 0.0; // Equivalent to std::numeric_limits<double>::max() in C++
                int local_bestCity = -1;
                int local_bestInsertionIndex = -1;
                double vn1_vk, vn2_vk, vn1_vn2;
                double insertionCost;

                for (int vn = 0; vn < tourSize; vn++) {
                    int vn_1 = tour[vn];
                    int vn_2 = tour[(vn + 1) % tourSize];  // Circular tour
                    vn1_vk = distanceMatrix[vn_1][vk];
                    vn2_vk = distanceMatrix[vn_2][vk];
                    vn1_vn2 = distanceMatrix[vn_1][vn_2];
                    insertionCost = vn1_vk + vn2_vk - vn1_vn2;

                    if (insertionCost < local_minCost) {
                        local_minCost = insertionCost;
                        local_bestCity = vk;
                        local_bestInsertionIndex = vn;
                    }
                }

                if (local_minCost < minCost) {
                    minCost = local_minCost;
                    bestCity = local_bestCity;
                    bestInsertionIndex = local_bestInsertionIndex;
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

int main() {
    struct timespec start_time, end_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);

    int numThreads = 12;
    omp_set_num_threads(numThreads);

    char const *fileName = "9_coords.coord";
    //char const *fileName = "4096_coords.coord";
    int numCoordinates = readNumOfCoords(fileName);
    double **coords = readCoords(fileName, numCoordinates);

    // Allocate memory for the distance matrix
    double **distanceMatrix = (double **)malloc(numCoordinates * sizeof(double *));
    for (int i = 0; i < numCoordinates; i++) {
        distanceMatrix[i] = (double *)malloc(numCoordinates * sizeof(double));
    }

    createDistanceMatrix(coords, numCoordinates, distanceMatrix);

    int tour[numCoordinates];
    cheapestInsertion(distanceMatrix, numCoordinates, tour);

    clock_gettime(CLOCK_MONOTONIC, &end_time);
    double execution_time = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_nsec - start_time.tv_nsec) / 1e9;
    printf("Execution time: %lf seconds\n", execution_time);

    // Print the tour
    printf("Tour order: ");
    for (int i = 0; i < numCoordinates; i++) {
        printf("%d", tour[i]);
        if (i < numCoordinates - 1) {
            printf(" -> ");
        }
    }
    printf(" -> 0\n");

    // Correct Solution for 9_coords
    printf("Tour order: 0 -> 2 -> 6 -> 1 -> 8 -> 7 -> 3 -> 5 -> 4 -> 0\n");

    // Free Memory
    for (int i = 0; i < numCoordinates; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);
    }

    free(coords);
    free(distanceMatrix);

    return 0;
}
