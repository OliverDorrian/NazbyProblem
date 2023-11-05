#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "coordReader.h"

double calculateDistance(double x1, double y1, double x2, double y2) {
    // Calculate the Euclidean distance between two points.
    return sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}

void createDistanceMatrix(double **passedCoords, int numCoordinates, double **distanceMatrix) {
    for (int i = 0; i < numCoordinates; ++i) {
        for (int j = i + 1; j < numCoordinates; ++j) {
            double distance = calculateDistance(passedCoords[i][0], passedCoords[i][1], passedCoords[j][0], passedCoords[j][1]);
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
        }
        // Set diagonal elements to 0
        distanceMatrix[i][i] = 0.000000;
    }
}

#include <float.h>
#include <stdbool.h>


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

        for (int vn = 0; vn < step; vn++) {
            for (int vk = 0; vk < numCoordinates; vk++) {
                if (!visited[vk]) {
                    double dist_vn_vk = distanceMatrix[tour[vn]][vk];
                    if (dist_vn_vk > maxDistance) {
                        maxDistance = dist_vn_vk;
                        farthestVertex = vk;
                    }
                }
            }
        }

        // Insert the farthestVertex into the tour at the position that minimizes tour length
        int positionToInsert = -1;
        double minInsertionCost = INFINITY;

        for (int vn = 0; vn < step; vn++) {
            int vn1 = tour[(vn + 1) % step];
            double insertionCost = distanceMatrix[tour[vn]][farthestVertex] + distanceMatrix[farthestVertex][vn1] - distanceMatrix[tour[vn]][vn1];
            if (insertionCost < minInsertionCost) {
                minInsertionCost = insertionCost;
                positionToInsert = vn;
            }
        }

        // Shift elements to make space for the new vertex
        for (int i = step; i > positionToInsert; i--) {
            tour[i] = tour[i - 1];
        }
        tour[positionToInsert + 1] = farthestVertex;
        visited[farthestVertex] = 1;
    }

    // Reverse the order of elements in the tour array to get the correct tour
    for (int i = 0; i < numCoordinates / 2; i++) {
        int temp = tour[i];
        tour[i] = tour[numCoordinates - i - 1];
        tour[numCoordinates - i - 1] = temp;
    }
}


int main(int argc, char *argv[]) {

    argv[1] = "9_coords.coord";
    char const *fileName = argv[1];

    argv[2] = "fiOut.dat";
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