/* Include guard */
#ifndef COMMON_H
#define COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define N_part 10
#define GAS
#define NVE

#define RNG gsl_rng_mt19937 /* RNG to use */
#define SEED 12345  /* Seed for the initialization RNG */

/* Simple power functions */
#define pow2(x) (x)*(x)
#define pow6(x) pow2(x)*pow2(x)*pow2(x)
#define pow12(x) pow6(x)*pow6(x)

/* L-J potential constants */
#define EPS4 0.0416 /* Epsilon * 4 */
#define SIGMA 3.4

#define M 39.948 /* Mass for the argon atom in amu */

/* Particles are represented by its x,y,z-cordinates and velocities.  */
typedef struct particle {
    double x, y, z; /* Position cordinates */
    double v_x, v_y, v_z; /* Velocity components */
} particle;

#ifdef GAS
double gas_box_size();
double box_size() {
    return gas_box_size();
}
#else
#ifdef LIQUID
double liquid_box_size();
double box_size() {
    return liquid_box_size();
}
#endif
#endif

#ifdef NVE
void init_NVE_velocities(particle *, double);
void init_velocities(particle * particles, double E_tot) {
    init_NVE_velocities(particles, E_tot);
}
#else
#ifdef NVT
void init_NVT_velocities(particle *, double);
void init_velocities(particle * particles, double T) {
    init_NVE_velocities(particles, E_tot);
}
#endif
#endif

#endif