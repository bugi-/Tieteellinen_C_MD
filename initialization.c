#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>

int init_cordinates(int *grid, int num_particles, int length) {
    int i;
    for (i = 0; i < length * length; i++) {
        grid[i] = 0;
    }
    gsl_rng * rand_gen = gsl_rng_alloc(gsl_rng_mt19937);
    /* Generate random positions for the particles */
    for (i = 0; i < num_particles; i++) {
        int x = gsl_rng_uniform_int (rand_gen, length);
        int y = gsl_rng_uniform_int (rand_gen, length);
        if (grid[y * length + x]) { /* Particle already in the selected position. */
            i--;
            continue;
        }
        grid[y * length + x] = 1;
    }
    return 0;
}

void print_grid(int *grid, int L) {
    int i, j;
    for (i = 0; i < L; i++) {
        for (j = 0; j < L; j++) {
            printf("%u ", grid[i * L + j]);
        }
        putchar('\n');
    }
    putchar('\n');
}

int main(void) {
    
    int N = 10;
    int L = 10;
    int grid[L*L];
    init_cordinates(grid, N, L);
    print_grid(grid, L);
}