#include "potentials/harmonic_potential.h"

HarmonicPotential::HarmonicPotential(double mass, double omega) : mass(mass), omega(omega) {
    k = mass * omega * omega;
}

double HarmonicPotential::V(const dVec& x) {
    double potential = 0.0;

    for (int ptcl_idx = 0; ptcl_idx < x.len(); ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            potential += x(ptcl_idx, axis) * x(ptcl_idx, axis);
        }
    }

    potential *= 0.5 * k;

    return potential;
}

dVec HarmonicPotential::gradV(const dVec& x) {
    return k * x;
}

double HarmonicPotential::laplacianV(const dVec& x) {
    return k;
}