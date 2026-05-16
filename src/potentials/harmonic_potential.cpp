#include "potentials/harmonic_potential.h"

HarmonicPotential::HarmonicPotential(double mass, double omega) : mass(mass), omega(omega) {
    k = mass * omega * omega;
}

double HarmonicPotential::V(const SingleVec& x) {
    double potential = 0.0;
    for (int axis = 0; axis < NDIM; ++axis) {
        potential += x[axis] * x[axis];
    }
    return 0.5 * k * potential;
}

SingleVec HarmonicPotential::gradV(const SingleVec& x) {
    SingleVec grad{};
    for (int axis = 0; axis < NDIM; ++axis) {
        grad[axis] = k * x[axis];
    }
    return grad;
}

double HarmonicPotential::laplacianV(const SingleVec& /* x */) {
    return k * NDIM;
}