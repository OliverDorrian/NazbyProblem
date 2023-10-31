#include <iostream>
#include "coordReader.c"
#include <limits> // Include the limits header for numeric_limits

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
        distanceMatrix[i][i] = 0;
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
        double minCost = std::numeric_limits<double>::max();

        for (int vk = 0; vk < numCoordinates; vk++) {
            if (unvisited[vk] != -1) {  // Check for unvisited cities
                double minInsertionCost = std::numeric_limits<double>::max();
                int bestInsertionPos = -1;

                for (int vn = 0; vn < tourSize; vn++) {
                    int vn_1 = tour[vn];
                    int vn_2 = tour[(vn + 1) % tourSize];  // Circular tour

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
    int numCoordinates = readNumOfCoords(fileName);
    double** coords = readCoords(fileName,numCoordinates);

    // Allocate memory for the distance matrix
    auto **distanceMatrix = (double **)malloc(numCoordinates * sizeof(double *));
    for (int i = 0; i < numCoordinates; i++) {
        distanceMatrix[i] = (double *)malloc(numCoordinates * sizeof(double ));
    }
    // generate distance Matrix
    createDistanceMatrix(coords, numCoordinates, distanceMatrix);

    int tour[numCoordinates];
    cheapestInsertion(distanceMatrix, numCoordinates, tour);

    double totalTourDistance = 0.0;
    for (int i = 0; i < numCoordinates - 1; i++) {
        totalTourDistance += distanceMatrix[tour[i]][tour[i + 1]];
    }
    std::cout << "Total Tour Distance: " << totalTourDistance << std::endl;

    // Print the tour
    std::cout << "Tour order: ";
    for (int i = 0; i < numCoordinates; i++) {
        std::cout << tour[i];
        if (i < numCoordinates - 1) {
            std::cout << " -> ";
        }
    }
    std::cout << " -> " << 0 << std::endl; // Add the starting point at the end to close the loop

    // Correct Solution
    std::cout << "FTour order:0 -> 2 -> 6 -> 1 -> 8 -> 7 -> 3 -> 5 -> 4 -> 0 " << std::endl;

    // Free Memory
    for (int i = 0; i < numCoordinates; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);
    }
    free(coords);
    free(distanceMatrix);

    return 0;
}
