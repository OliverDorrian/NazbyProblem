#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "coordReader.c"

double calculateDistance(double x1, double y1, double x2, double y2) {
    // Calculate the Euclidean distance between two points.
    return sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}

void createDistanceMatrix(double **passedCoords, int numCoordinates, double **distanceMatrix) {

    for (int i = 0; i < numCoordinates; ++i) {
        for (int j = i + 1; j < numCoordinates; ++j) {
            double distance = calculateDistance(passedCoords[i][0], passedCoords[i][1], passedCoords[j][0],
                                                passedCoords[j][1]);
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
        }
        // Set diagonal elements to 0
        distanceMatrix[i][i] = 0.000000;
    }
}

void cheapestInsertion(double** distanceMatrix, int numCoordinates, int* tour) {
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
        double minCost = 1.79769e+308;

        for (int vk = 0; vk < numCoordinates; vk++) {
            if (unvisited[vk] != -1) {  // Check for unvisited cities
                double minInsertionCost = 1.79769e+308;
                int bestInsertionPos = -1;

                for (int vn = 0; vn < tourSize; vn++) {
                    int vn_1 = tour[vn];
                    int vn_2 = tour[(vn + 1) % tourSize];

                    double insertionCost = distanceMatrix[vn_1][vk] + distanceMatrix[vn_2][vk] - distanceMatrix[vn_1][vn_2];

                    if (insertionCost < minInsertionCost) {
                        minInsertionCost = insertionCost;
                        bestInsertionPos = vn;
                    }
                }

                if (minInsertionCost < minCost) {
                    minCost = minInsertionCost;
                    bestCity = vk;
                    bestInsertionIndex = bestInsertionPos;
                }
            }
        }

        tourSize++;  // Increment tour size
        int nextTourIndex = (bestInsertionIndex + 1) % tourSize;  // Circular tour
        for (int vn = tourSize - 1; vn > nextTourIndex; vn--) {
            tour[vn] = tour[vn - 1];
        }
        tour[nextTourIndex] = bestCity;
        unvisited[bestCity] = -1;  // Mark the chosen city as visited
    }
}

int main() {
    char const *fileName = "9_coords.coord";
    //char const *fileName = "4096_coords.coord";

    int numCoordinates = readNumOfCoords(fileName);
    double** coords = readCoords(fileName,numCoordinates);

    // Allocate memory for the distance matrix
    double **distanceMatrix = (double **)malloc(numCoordinates * sizeof(double *));
    for (int i = 0; i < numCoordinates; i++) {
        distanceMatrix[i] = (double *)malloc(numCoordinates * sizeof(double ));
    }

    clock_t start_time = clock();

    createDistanceMatrix(coords, numCoordinates, distanceMatrix);

    int tour[numCoordinates];
    cheapestInsertion(distanceMatrix, numCoordinates, tour);


    clock_t end_time = clock();
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution time: %.6f seconds\n", execution_time);

    // Print the tour
    printf("Tour order: ");
    for (int i = 0; i < numCoordinates; i++) {
        printf("%d", tour[i]);
        if (i < numCoordinates - 1) {
            printf(" -> ");
        }
    }
    printf(" -> 0\n"); // Add the starting point at the end to close the loop

    // Correct Solution
    printf("FTour order:0 -> 2 -> 6 -> 1 -> 8 -> 7 -> 3 -> 5 -> 4 -> 0\n");

    // Free Memory
    for (int i = 0; i < numCoordinates; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);
    }
    free(coords);
    free(distanceMatrix);



    return 0;
}
