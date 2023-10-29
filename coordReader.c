#include <cstdio>
#include <cstdlib>
#include <cmath>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wwritable-strings"
#pragma clang diagnostic pop

int readNumOfCoords(const char *fileName);
double **readCoords(const char *filename, int numOfCoords);

/*Gets the number of the coordinates in the file. Returns as a single integer*/
int readNumOfCoords(const char *filename){
    FILE *file = fopen(filename, "r");
    int numOfCoords = 0;

    if(file == nullptr){
        return -1;
    }

    char line[100];

    while(fgets(line, sizeof(line), file) != nullptr){
        numOfCoords++;
    }

    return numOfCoords;
}

/*Gets the data from the file. Returns as an array of doubles, Ignores the first integer*/
double **readCoords(const char *filename, int numOfCoords){
    FILE *file = fopen(filename,"r");
    int i;

    char line[100];

    if(file == nullptr) {
        printf("Unable to open file: %s", filename);
        return nullptr;
    }

    double **coords = (double **)malloc(numOfCoords * sizeof(double *));

    for(i = 0; i < numOfCoords; i++){
        coords[i] = (double *) malloc(2 * sizeof(double));
        if (coords[i] == nullptr){
            perror("Memory Allocation Failed");
        }
    }

    int lineNum = 0;
    while(fgets(line, sizeof(line), file) != nullptr){
        double x, y;
        if (sscanf(line, "%lf,%lf", &x, &y) == 2){
            coords[lineNum][0] = x;
            coords[lineNum][1] = y;
            lineNum++;
        }
    }

    return coords;
}

void *writeTourToFile(int *tour, int tourLength, char *filename){

    FILE *file = fopen(filename, "w");
    int i;

    if(file == nullptr){
        printf("Unable to open file: %s", filename);
        return nullptr;
    }

    fprintf(file, "%d \n", tourLength);

    printf("Writing output data\n");
    for(i=0; i < tourLength; i++) {
        fprintf(file, "%d ", tour[i]);
    }
    return nullptr;
}
