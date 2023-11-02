#include <iostream>
#include <limits>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <numeric>
#include <vector>
#include <list>
#include <algorithm>

#include "coordReader.c"

double calculateDistance(double x1, double y1, double x2, double y2) {
    // Calculate the Euclidean distance between two points.
    return sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}

void checkTreads() {
    #pragma omp parallel
    {
        int threadID = omp_get_thread_num();
    #pragma omp critical
        std::cout << "Hello from thread " << threadID << std::endl;
    }
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

void printDistanceMatrix(double **distanceMatrix, int numCoordinates) {

    for (int i = 0; i < numCoordinates; i++) {
        for (int j = 0; j < numCoordinates; j++) {
            std::cout << distanceMatrix[i][j] << " ";
        }
        std::cout << "\n";
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
        double minCost = std::numeric_limits<double>::max();

        #pragma omp parallel for
        for (int vk = 0; vk < numCoordinates; vk++) {
            if (unvisited[vk] != -1) {  // Check for unvisited cities
                double local_minCost = std::numeric_limits<double>::max();
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

                if (local_minCost < minCost) {
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
    auto start_time = std::chrono::high_resolution_clock::now();

    int numThreads = 12;
    omp_set_num_threads(numThreads);

    //char const *fileName = "9_coords.coord";
    char const *fileName = "4096_coords.coord";
    int numCoordinates = readNumOfCoords(fileName);
    double **coords = readCoords(fileName, numCoordinates);


    // Allocate memory for the distance matrix
    auto **distanceMatrix = (double **) malloc(numCoordinates * sizeof(double *));
    for (int i = 0; i < numCoordinates; i++) {
        distanceMatrix[i] = (double *) malloc(numCoordinates * sizeof(double));
    }

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
    std::cout << " -> " << 0 << std::endl;

    // Correct Solution for 9_coords
    std::cout << "Tour order: 0 -> 2 -> 6 -> 1 -> 8 -> 7 -> 3 -> 5 -> 4 -> 0 " << std::endl;

    // Free Memory
    for (int i = 0; i < numCoordinates; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);
    }

    free(coords);
    free(distanceMatrix);

    return 0;
}