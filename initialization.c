#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define RNG gsl_rng_mt19937 /* RNG to use */
#define SEED 12345  /* Seed for the initialization RNG */
#define DELTA_PART 0.5

gsl_rng *rand_gen;

/* Particles are represented by its x,y-cordinates.  */
typedef struct particle {
    double x, y, z; /* Position cordinates */
    double v_x, v_y, v_z; /* Velocity components */
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
    for (int i = 0; i < num_particles; i++) {
        /* Generate x, y and z for a particle */
        double rand1 = length * gsl_rng_uniform(rand_gen);
        double rand2 = length * gsl_rng_uniform(rand_gen);
        double rand3 = length * gsl_rng_uniform(rand_gen);
        /* Check if a particle with these coordinates already exists */
        int found = 0;
        for (int j = 0; j < i; j++) {
            if (particles[j].x - rand1 < DELTA_PART && particles[j].y - rand2 < DELTA_PART && particles[j].z - rand3 < DELTA_PART) {
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

void init_NVE_velocities(particle * particles, int N) {
    for (int i = 0; i < N; i++) {
        particles[i].v_x = 1;
    }
}

void init_NVT_velocities() {
    
}

double potential_energy(int *grid) {
    return 0.0;
}

int main(void) {
    rand_gen = gsl_rng_alloc(RNG);
    gsl_rng_set(rand_gen, SEED);
    int N = 10;
    particle *particles = (particle *)calloc(N, sizeof(particle));
    init_gas_particles(particles, N);
    for (int i = 0; i < N; i++) {
        printf("%f %f %f\n", particles[i].x, particles[i].y, particles[i].z);
    }
    return 0;
    
    gsl_rng_free(rand_gen);
}