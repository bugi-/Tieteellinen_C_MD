/* Include guard */
#ifndef COMMON_H
#define COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define N_part 10 /* Number of particles */
#define GAS /* 'GAS' or 'LIQUID' */
#define NVE /* 'NVE' or 'NVT' */

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

typedef struct {
    double x, y, z;
} vector;

void remove_CM_velocity(particle *);

double box_size();

void init_velocities(particle * , double);

#endif