#ifndef COMMON_H
#include "common.h"
#endif

void check_boundaries(particle *particles) {
    double L = box_size();
    for (int i = 0; i < N_part; i++) {
        if (particles[i].x > L) particles[i].x -= L;
        if (particles[i].y > L) particles[i].y -= L;
        if (particles[i].z > L) particles[i].z -= L;
        if (particles[i].x < 0) particles[i].x += L;
        if (particles[i].y < 0) particles[i].y += L;
        if (particles[i].z < 0) particles[i].z += L;
    }
}
