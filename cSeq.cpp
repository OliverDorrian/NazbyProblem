#include <iostream>
#include "coordReader.c"
#include <limits>
#include <chrono>
#include <iomanip>
#include <numeric>

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

void printDistanceMatrix(double **distanceMatrix, int numCoordinates) {

    for (int i = 0; i < numCoordinates; i++) {
        for (int j = 0; j < numCoordinates; j++) {
            std::cout << std::fixed << std::setprecision(2) << std::setw(6) << distanceMatrix[i][j] << " ";
        }
        std::cout << "\n";
    }
}

void cheapestInsertion(double** distanceMatrix, int numCoordinates, int* tour) {
    int unvisited[numCoordinates];
    std::iota(unvisited, unvisited + numCoordinates, 0);

    int tourSize = 1;
    int startCity = 0;
    tour[0] = startCity;
    unvisited[startCity] = -1;  // Mark the starting city as visited

    for (int i = 1; i < numCoordinates; i++) {
        std::cout << i << std::endl;
        int bestCity = -1;
        int bestInsertionIndex = -1;
        double minCost = std::numeric_limits<double>::max();

        for (int vk = 0; vk < numCoordinates; vk++) {
            if (unvisited[vk] != -1) {  // Check for unvisited cities
                double minInsertionCost = std::numeric_limits<double>::max();
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
    auto **distanceMatrix = (double **)malloc(numCoordinates * sizeof(double *));
    for (int i = 0; i < numCoordinates; i++) {
        distanceMatrix[i] = (double *)malloc(numCoordinates * sizeof(double ));
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    createDistanceMatrix(coords, numCoordinates, distanceMatrix);

    int tour[numCoordinates];
    cheapestInsertion(distanceMatrix, numCoordinates, tour);


    auto end_time = std::chrono::high_resolution_clock::now();
    auto execution_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    double seconds = static_cast<double>(execution_time.count()) / 1000000.0;
    std::cout << "Execution time: " << seconds << " seconds" << std::endl;

    end_time = std::chrono::high_resolution_clock::now();
    execution_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "Execution time: " << execution_time.count() << " miroseconds" << std::endl;

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
