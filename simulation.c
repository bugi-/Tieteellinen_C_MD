#include "common.h"

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

double kinetic_energy(particle *particles) {
    double E_k = 0.0;
    for (int i = 0; i < N_part; i++) {
        double v = pow2(particles[i].v_x) + pow2(particles[i].v_y) + pow2(particles[i].v_z);
        E_k += v;
    }
    E_k *= 0.5 * 103.6427 * M;
    return E_k;
}

void remove_angular_moment(particle *particles) {
    
}

vector* calc_force() {
    return NULL;
}

void timestep(particle *particles) {
    
}