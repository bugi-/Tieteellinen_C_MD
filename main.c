#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

#define NUM_ITER 100000000

int main(void) {
int i;
    long N_accepted = 0;
    gsl_rng *rand_gen;
        rand_gen = gsl_rng_alloc(gsl_rng_mt19937);

    for (i = 0; i < NUM_ITER; i++) {
        double x = gsl_rng_uniform(rand_gen);
        double y = gsl_rng_uniform(rand_gen);
        if (x*x + y*y <= 1) {
            N_accepted++;
        }
    }

    printf("Estimate for pi with N=%d: %f\n", 4.0 * N_accepted / NUM_ITER, NUM_ITER);

    return 0;
}
