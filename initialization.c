#ifndef COMMON_H
#include "common.h"
#endif
#include <gsl/gsl_randist.h>

#define EPS_PART 0.5 /* Epsilon for closeness check of particles*/

gsl_rng *rand_gen;

/* Define length of a side in the simulation grid in Å. Sizes are the ones asked to be used. */
/* Cast to int to make things easier */
double gas_box_size(){
    return 33.38 * cbrt(N_part);
}
double liquid_box_size(){
    return 3.62 * cbrt(N_part);
}

/* Calculates E_pot for the system */
double potential_energy(particle *particles) {
    double tot_U = 0.0;
    double half_L = 0.5 * box_size();
    for (int i = 0; i < N_part; i++) {
        for (int j = 0; j < i; j++) {
            double dx = abs(particles[i].x - particles[j].x);
            double dy = abs(particles[i].y - particles[j].y);
            double dz = abs(particles[i].z - particles[j].z);
            /* Boundary conditions */
            if (dx > half_L) dx -= half_L;
            if (dy > half_L) dy -= half_L;
            if (dz > half_L) dz -= half_L;
            double r_ij = sqrt(pow2(dx) + pow2(dy) + pow2(dz));
            tot_U += EPS4 * (pow12(SIGMA/r_ij) - pow6(SIGMA/r_ij));
        }
    }
    return tot_U;
}

void init_particles(particle *particles) {
    double length = box_size();
    for (int i = 0; i < N_part; i++) {
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

//void init_gas_particles(particle *particles) {
//    initialize_particles(particles);
//}
//
//void init_liquid_particles(particle *particles) {
//    initialize_particles(particles);
//}

/* Functions for initializing velocities for the particles in both cases. Particles have to be generated for these to work. */
void init_NVE_velocities(particle * particles, double E_tot) {
    double gauss_nums[N_part];
    double sum_gauss = 0;
    double v_components[3];
    for (int i = 0; i < N_part; i++) {
        gauss_nums[i] = gsl_ran_ugaussian(rand_gen);
        sum_gauss += gauss_nums[i];
    }
    double E_kin = E_tot - potential_energy(particles);
    E_kin *= 2;
    E_kin /= M;
    for (int i = 0; i < N_part; i++) {
        gauss_nums[i] /= sum_gauss;
        gauss_nums[i] = sqrt(gauss_nums[i]); /* Use the same array for speeds to save space */
        v_components[1] = gsl_rng_uniform(rand_gen);
        v_components[2] = gsl_rng_uniform(rand_gen);
        v_components[3] = gsl_rng_uniform(rand_gen);
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

void init_NVT_velocities(particle * particles, double T) {
    for (int i = 0; i < N_part; i++) {
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
    particle *particles = (particle *)calloc(N_part, sizeof(particle));
    init_particles(particles);
    init_velocities(particles, .0259);
    for (int i = 0; i < N_part; i++) {
        printf("%f %f %f %f %f %f\n", particles[i].x, particles[i].y, particles[i].z, particles[i].v_x, particles[i].v_y, particles[i].v_z);
    }
    gsl_rng_free(rand_gen);
    return 0;
}