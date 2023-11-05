#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <float.h>
#include "coordReader.h"

double calculateDistance(double x1, double y1, double x2, double y2) {
    // Calculate the Euclidean distance between two points.
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

void farthestInsertion(double **distanceMatrix, int numCoordinates, int *tour) {
    int visited[numCoordinates];
    for (int i = 0; i < numCoordinates; i++) {
        visited[i] = 0;
        tour[i] = -1;
    }

    // Start with the first vertex as the initial tour
    int currentVertex = 0;
    tour[0] = currentVertex;
    visited[currentVertex] = 1;

    for (int step = 1; step < numCoordinates; step++) {
        int farthestVertex = -1;
        double maxDistance = -1.0;

        #pragma omp parallel for
        for (int vk = 0; vk < numCoordinates; vk++) {
            if (!visited[vk]) {
                double localMaxDistance = -1.0;
                int localFarthestVertex = -1;
                for (int vn = 0; vn < step; vn++) {
                    double dist_vn_vk = distanceMatrix[tour[vn]][vk];
                    if (dist_vn_vk > localMaxDistance) {
                        localMaxDistance = dist_vn_vk;
                        localFarthestVertex = vk;
                    }
                }
                #pragma omp critical
                {
                    if (localMaxDistance > maxDistance) {
                        maxDistance = localMaxDistance;
                        farthestVertex = localFarthestVertex;
                    }
                }
            }
        }

        // Insert the farthestVertex into the tour at the position that minimizes tour length
        int positionToInsert = -1;
        double minInsertionCost = DBL_MAX;

        #pragma omp parallel for
        for (int vn = 0; vn < step; vn++) {
            int vn1 = tour[(vn + 1) % step];
            double insertionCost = distanceMatrix[tour[vn]][farthestVertex] + distanceMatrix[farthestVertex][vn1] - distanceMatrix[tour[vn]][vn1];
            #pragma omp critical
            {
                if (insertionCost < minInsertionCost) {
                    minInsertionCost = insertionCost;
                    positionToInsert = vn;
                }
            }
        }

        // Shift elements to make space for the new vertex

        for (int i = step; i > positionToInsert; i--) {
            tour[i] = tour[i - 1];
        }
        tour[positionToInsert + 1] = farthestVertex;
        visited[farthestVertex] = 1;
    }
}


int main(int argc, char *argv[]) {
    int numThreads = 12;
    omp_set_num_threads(numThreads);
    argv[1] = "9_coords.coord";
    //argv[1] = "4096_coords.coord";
    char const *fileName = argv[1];

    argv[2] = "icompOut.dat";
    //char const *fileName = "4096_coords.coord";
    //char const *fileName = "9_coords.coord";

    int numCoordinates = readNumOfCoords(fileName);
    double** coords = readCoords(fileName,numCoordinates);

    clock_t start_time = clock();

    double **distanceMatrix = (double **)malloc(numCoordinates * sizeof(double *));
    for (int i = 0; i < numCoordinates; i++) {
        distanceMatrix[i] = (double *)malloc(numCoordinates * sizeof(double));
    }

    // Generate distance Matrix
    createDistanceMatrix(coords, numCoordinates, distanceMatrix);

    int tour[numCoordinates];
    farthestInsertion(distanceMatrix, numCoordinates, tour);

    clock_t end_time = clock();
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution time: %.6f seconds\n", execution_time);

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