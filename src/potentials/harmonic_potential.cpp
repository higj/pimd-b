#include "potentials/harmonic_potential.h"

HarmonicPotential::HarmonicPotential(double mass, double omega) : mass(mass), omega(omega) {
    k = mass * omega * omega;
}

double HarmonicPotential::V(const Vec& x) {
    double potential = 0.0;
    for (int axis = 0; axis < NDIM; ++axis) {
        potential += x[axis] * x[axis];
    }
    return 0.5 * k * potential;
}

Vec HarmonicPotential::gradV(const Vec& x) {
    Vec grad{};
    for (int axis = 0; axis < NDIM; ++axis) {
        grad[axis] = k * x[axis];
    }
    return grad;
}

double HarmonicPotential::laplacianV(const Vec& /* x */) {
    return k * NDIM;
}