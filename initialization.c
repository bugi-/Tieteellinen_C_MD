#include "common.h"
#include <gsl/gsl_randist.h>

#define EPS_PART 0.5 /* Epsilon for closeness check of particles*/

gsl_rng *rand_gen;

/* Define length of a side in the simulation grid in Å. Sizes are the ones asked to be used. */
/* Cast to int to make things easier */
int gas_box_size(int N){
    return (int)(33.38 * cbrt(N));
}
int liquid_box_size(int N){
    return (int)(3.62 * cbrt(N));
}

/* Calculates E_pot for the system */
double potential_energy(particle *particles, int N) {
    double tot_U = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            double r_ij = sqrt(pow2(abs(particles[i].x - particles[j].x))
            + pow2(abs(particles[i].y - particles[j].y))
            + pow2(abs(particles[i].z - particles[j].z)));
            tot_U += EPS4 * (pow12(SIGMA/r_ij) - pow6(SIGMA/r_ij));
        }
    }
    return tot_U;
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
            if (particles[j].x - rand1 < EPS_PART && particles[j].y - rand2 < EPS_PART && particles[j].z - rand3 < EPS_PART) {
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

/* Functions for initializing velocities for the particles in both cases. Particles have to be generated for these to work. */
void init_NVE_velocities(particle * particles, int N, double E_tot) {
    double gauss_nums[N];
    double sum_gauss = 0;
    double v_components[3];
    for (int i = 0; i < N; i++) {
        gauss_nums[i] = gsl_ran_ugaussian(rand_gen);
        sum_gauss += gauss_nums[i];
    }
    double E_kin = E_tot - potential_energy(particles, N);
    E_kin *= 2;
    E_kin /= M;
    for (int i = 0; i < N; i++) {
        gauss_nums[i] /= sum_gauss;
        gauss_nums[i] = sqrt(gauss_nums[i]); /* Use the same array for speeds to save space */
        v_components[1] = gsl_ran_ugaussian(rand_gen);
        v_components[2] = gsl_ran_ugaussian(rand_gen);
        v_components[3] = gsl_ran_ugaussian(rand_gen);
        double sum_v_comp = pow2(v_components[1]) + pow2(v_components[2]) + pow2(v_components[3]);
        sum_v_comp = sqrt(sum_v_comp);
        v_components[1] /= sum_v_comp;
        v_components[2] /= sum_v_comp;
        v_components[3] /= sum_v_comp;
        particles[i].v_x = v_components[1] * gauss_nums[i];
        particles[i].v_y = v_components[2] * gauss_nums[i];
        particles[i].v_z = v_components[3] * gauss_nums[i];
    }
    
}

void init_NVT_velocities(particle * particles, int N, double T) {
    for (int i = 0; i < N; i++) {
        double rand1 = gsl_ran_ugaussian(rand_gen);
        double rand2 = gsl_ran_ugaussian(rand_gen);
        double rand3 = gsl_ran_ugaussian(rand_gen);
        double scale_factor = sqrt(T / M);
        particles[i].v_x = rand1 * scale_factor;
        particles[i].v_y = rand2 * scale_factor;
        particles[i].v_z = rand3 * scale_factor;
    }
}

int main(void) {
    rand_gen = gsl_rng_alloc(RNG);
    gsl_rng_set(rand_gen, SEED);
    int N = 10;
    particle *particles = (particle *)calloc(N, sizeof(particle));
    init_gas_particles(particles, N);
//    for (int i = 0; i < N; i++) {
//        printf("%f %f %f\n", particles[i].x, particles[i].y, particles[i].z);
//    }
    gsl_rng_free(rand_gen);
    return 0;
}