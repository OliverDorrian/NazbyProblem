#include <iostream>
#include <limits>
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
        distanceMatrix[i][i] = 0;
    }
}

void farthestInsertion(double** distanceMatrix, int numCoordinates, int* tour) {
    bool visited[numCoordinates];
    for (int i = 0; i < numCoordinates; i++) {
        visited[i] = false;
        tour[i] = -1;
    }

    // Start with the first vertex as the initial tour
    int currentVertex = 0;
    tour[0] = currentVertex;
    visited[currentVertex] = true;

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
        double minInsertionCost = std::numeric_limits<double>::max();

        for (int vn = 0; vn < step; vn++) {
            int vn1 = tour[(vn + 1) % step];
            double insertionCost = distanceMatrix[tour[vn]][farthestVertex] + distanceMatrix[farthestVertex][vn1] - distanceMatrix[tour[vn]][vn1];
            if (insertionCost < minInsertionCost) {
                minInsertionCost = insertionCost;
                positionToInsert = vn;
            }
        }

        // Insert the farthestVertex at the appropriate position in the tour
        for (int i = step; i > positionToInsert; i--) {
            tour[i] = tour[i - 1];
        }
        tour[positionToInsert + 1] = farthestVertex;
        visited[farthestVertex] = true;
    }
    // Connect the last vertex to the starting vertex to form a closed tour
    tour[numCoordinates] = tour[0];
}



int main () {
    char const *fileName = "9_coords.coord";
    int numCoordinates = readNumOfCoords(fileName);
    double** coords = readCoords(fileName,numCoordinates);

    auto **distanceMatrix = (double **)malloc(numCoordinates * sizeof(double *));
    for (int i = 0; i < numCoordinates; i++) {
        distanceMatrix[i] = (double *)malloc(numCoordinates * sizeof(double ));
    }

    // generate distance Matrix
    createDistanceMatrix(coords, numCoordinates, distanceMatrix);

    int tour[numCoordinates];
    farthestInsertion(distanceMatrix, numCoordinates, tour);

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

    // Correct Solution String
    std::cout << "FTour order:0 -> 4 -> 2 -> 5 -> 3 -> 7 -> 8 -> 1 -> 6 -> 0" << std::endl;

    // Free Memory
    for (int i = 0; i < numCoordinates; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);
    }
    free(coords);
    free(distanceMatrix);

}