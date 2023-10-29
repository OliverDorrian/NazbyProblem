#include <iostream>
#include "coordReader.c"

double calculateDistance(double x1, double y1, double x2, double y2) {
    // Calculate the Euclidean distance between two points.
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

int main() {
    char const *fileName = "9_coords.coord";

    int numCoordinates = readNumOfCoords(fileName);
    double** coords;
    coords = readCoords(fileName,numCoordinates);

    int* visited = (int*)malloc(numCoordinates * sizeof(int));
    int tour[numCoordinates];

    for (int i = 0; i < numCoordinates; i++) {
        visited[i] = 0;  // Mark all coordinates as unvisited
    }

    int current = 0;  // Start from the first coordinate
    tour[0] = current;
    visited[current] = 1;  // Mark the first coordinate as visited

    for (int i = 1; i < numCoordinates; i++) {
        int nearest = -1;
        double minDistance = -1;

        for (int j = 0; j < numCoordinates; j++) {
            if (!visited[j]) {
                double distance = calculateDistance(coords[current][0], coords[current][1], coords[j][0], coords[j][1]);
                if (nearest == -1 || distance < minDistance) {
                    nearest = j;
                    minDistance = distance;
                }
            }
        }

        tour[i] = nearest;
        visited[nearest] = 1;
        current = nearest;
    }

    return 0;
}
