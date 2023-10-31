#include <iostream>
#include "coordReader.c"

double calculateDistance(double x1, double y1, double x2, double y2) {
    // Calculate the Euclidean distance between two points.
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
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

void cheapestInsertion(double** distanceMatrix, int numCoordinates, int* tour) {
    int unvisited[numCoordinates];
    for (int i = 0; i < numCoordinates; i++) {
        unvisited[i] = i;
    }

    int tourSize = 0;
    int startCity = 0;
    tour[0] = startCity;
    unvisited[startCity] = -1;

    while (tourSize < numCoordinates) {
        int bestCity = -1;
        int bestInsertionIndex = -1;
        double minCost = std::numeric_limits<double>::max();

        for (int i = 0; i <= tourSize; i++) {
            for (int j = 0; j < numCoordinates; j++) {
                if (unvisited[j] != -1) {
                    double cost = distanceMatrix[tour[i]][j];
                    if (cost < minCost) {
                        minCost = cost;
                        bestCity = j;
                        bestInsertionIndex = i;
                    }
                }
            }
        }

        tourSize++;
        tour[tourSize] = bestCity;
        unvisited[bestCity] = -1;

        // Shift the tour to insert the bestCity
        for (int i = tourSize; i > bestInsertionIndex; i--) {
            tour[i] = tour[i - 1];
        }
        tour[bestInsertionIndex] = bestCity;
    }
}

int main() {
    char const *fileName = "9_coords.coord";
    int numCoordinates = readNumOfCoords(fileName);
    double** coords = readCoords(fileName,numCoordinates);

    // Allocate memory for the distance matrix
    auto **distanceMatrix = (double **)malloc(numCoordinates * sizeof(double *));
    for (int i = 0; i < numCoordinates; i++) {
        distanceMatrix[i] = (double *)malloc(numCoordinates * sizeof(double ));
    }
    // generate distance Matrix
    createDistanceMatrix(coords, numCoordinates, distanceMatrix);

    int* visited = (int*)malloc(numCoordinates * sizeof(int));
    int tour[numCoordinates];


    // insert additions here
    cheapestInsertion(distanceMatrix, numCoordinates, tour);

    // Print the tour
    std::cout << "Tour order: ";
    for (int i = 0; i < numCoordinates; i++) {
        std::cout << tour[i];
        if (i < numCoordinates - 1) {
            std::cout << " -> ";
        }
    }
    std::cout << " -> " << tour[0] << std::endl; // Add the starting point at the end to close the loop

    // Correct Solution
    std::cout << "CTour order:0 -> 2 -> 6 -> 1 -> 8 -> 7 -> 3 -> 5 -> 4 -> 0 " << std::endl;

    // Free Memory
    for (int i = 0; i < numCoordinates; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);
    }
    free(coords);
    free(distanceMatrix);

    return 0;
}
