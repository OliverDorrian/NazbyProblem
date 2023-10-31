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


int main() {
    char const *fileName = "9_coords.coord";
    int numCoordinates = readNumOfCoords(fileName);
    double** coords = readCoords(fileName,numCoordinates);

    // Allocate memory for the distance matrix
    auto **distanceMatrix = (double **)malloc(numCoordinates * sizeof(double *));
    for (int i = 0; i < numCoordinates; i++) {
        distanceMatrix[i] = (double *)malloc(numCoordinates * sizeof(double ));
    }

    createDistanceMatrix(coords, numCoordinates, distanceMatrix);


    int* visited = (int*)malloc(numCoordinates * sizeof(int));
    int tour[numCoordinates];

    for (int i = 0; i < numCoordinates; i++) {
        visited[i] = 0;  // Mark all coordinates as unvisited
    }

    int current = 0;  // Start from the first coordinate
    tour[0] = current;
    visited[current] = 0;  // Mark the first coordinate as visited

    for (int i = 0; i < numCoordinates; i++) {
        int nearest = -1;
        double minDistance = -1;

        for (int j = 0; j < numCoordinates; j++) {
            if (!visited[j]){
                for (int k = 0; k < ; ++k) {
                    
                }
            }
            (nearest == -1 || distanceMatrix[current][j] < minDistance)) {
                nearest = j;
                minDistance = distanceMatrix[current][j];
            }
        }

        tour[i] = nearest;
        visited[nearest] = 1;
        current = nearest;
    }

    // Print the tour
    std::cout << "Tour order: ";
    for (int i = 0; i < numCoordinates; i++) {
        std::cout << tour[i];
        if (i < numCoordinates - 1) {
            std::cout << " -> ";
        }
    }

    std::cout << " -> " << tour[0] << std::endl; // Add the starting point at the end to close the loop


    std::cout << "CTour order:0 -> 2 -> 6 -> 1 -> 8 -> 7 -> 3 -> 5 -> 4 -> 0 " << std::endl;

    for (int i = 0; i < numCoordinates; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);
    }
    free(coords);
    free(distanceMatrix);

    return 0;
}
