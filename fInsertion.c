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
        distanceMatrix[i][i] = 0;
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
}

int main() {

    char const *fileName = "4096_coords.coord";

    int numCoordinates = readNumOfCoords(fileName);
    double** coords = readCoords(fileName,numCoordinates);

    double **distanceMatrix = (double **)malloc(numCoordinates * sizeof(double *));
    for (int i = 0; i < numCoordinates; i++) {
        distanceMatrix[i] = (double *)malloc(numCoordinates * sizeof(double));
    }

    // Generate distance Matrix
    createDistanceMatrix(coords, numCoordinates, distanceMatrix);

    int tour[numCoordinates];
    farthestInsertion(distanceMatrix, numCoordinates, tour);

    double totalTourDistance = 0.0;
    for (int i = 0; i < numCoordinates - 1; i++) {
        totalTourDistance += distanceMatrix[tour[i]][tour[i + 1]];
    }
    printf("Total Tour Distance: %f\n", totalTourDistance);

    // Print the tour
    printf("Tour order: ");
    for (int i = 0; i < numCoordinates; i++) {
        printf("%d", tour[i]);
        if (i < numCoordinates - 1) {
            printf(" -> ");
        }
    }
    printf(" -> 0\n"); // Add the starting point at the end to close the loop

    // Correct Solution String
    printf("FTour order:0 -> 4 -> 2 -> 5 -> 3 -> 7 -> 8 -> 1 -> 6 -> 0\n");

    // Free Memory
    for (int i = 0; i < numCoordinates; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);
    }
    free(coords);
    free(distanceMatrix);
    return 0;
}