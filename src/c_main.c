#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main() {
    double coords[2][3];
    unsigned int i, j;
    for (i=0; i<2; i++) {
        for (j=0; j<3; j++) {
            coords[i][j] = 3 * i + j;
            printf("coords() %i %i %f\n", i + 1, j + 1, coords[i][j]);
        }
    }

    print_coords(coords); 

    return 0;
}
