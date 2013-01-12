#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

/* Particles are represented by its x,y-cordinates.  */
typedef struct partic {
    int x, y, z;
    int v_x, v_y, v_z;
} particle;

/* Define length of a side in the simulation grid in Å. Sizes are the ones asked to be used. */
/* Cast to int to make things easier */
int gas_box_size(int N){
    return (int)(33.38 * cbrt(N));
}
int liquid_box_size(int N){
    return (int)(3.62 * cbrt(N));
}

void init_particles(particle *particles, int num_particles, int length) {
    gsl_rng *rand_gen = gsl_rng_alloc(gsl_rng_mt19937);
    for (int i = 0; i < num_particles; i++) {
        /* Generate x and y for a particle */
        int rand1 = gsl_rng_uniform_int (rand_gen, length);
        int rand2 = gsl_rng_uniform_int (rand_gen, length);
        int rand3 = gsl_rng_uniform_int (rand_gen, length);
        /* Check if a particle with these coordinates already exists */
        int found = 0;
        for (int j = 0; j < i; j++) {
            if (particles[j].x == rand1 && particles[i].y == rand2 && particles[i].z == rand3) {
                i--;
                found = 1;
                break;
            }
        }
        if (!found) {
            particle new;
            new.x = rand1;
            new.y = rand2;
            new.z = rand3;
            particles[i] = new;
        }
    }
}

void init_gas_particles(particle *particles, int num_particles) {
    init_particles(particles, num_particles, gas_box_size(num_particles));
}

void init_liquid_particles(particle *particles, int num_particles) {
    init_particles(particles, num_particles, liquid_box_size(num_particles));
}

double potential_energy(int *grid) {
    return 0.0;
}

int main(void) {
    int i;
    int N = 10;
    particle *particles = (particle *)calloc(N, sizeof(particle));
    init_gas_particles(particles, N);
    for (i = 0; i < N; i++) {
        printf("%d %d\n", particles[i].x, particles[i].y);
    }
    return 0;
}