#include "potentials/double_well_potential.h"

DoubleWellPotential::DoubleWellPotential(double mass, double strength, double loc)
    : mass(mass), strength(strength), loc(loc) {}

double DoubleWellPotential::V(const SingleVec& x) {
    double potential = 0.0;
    const double loc2 = loc * loc;

    for (int axis = 0; axis < NDIM; ++axis) {
        potential += (x[axis] * x[axis] - loc2) * (x[axis] * x[axis] - loc2);
    }

    return potential * mass * strength;
}

SingleVec DoubleWellPotential::gradV(const SingleVec& x) {
    SingleVec grad{};
    const double loc2 = loc * loc;

    double prefactor = 0.0;
    for (int axis = 0; axis < NDIM; ++axis) {
        prefactor += x[axis] * x[axis];
    }
    prefactor = 4 * mass * strength * (prefactor - loc2);

    for (int axis = 0; axis < NDIM; ++axis) {
        grad[axis] = prefactor * x[axis];
    }

    return grad;
}

double DoubleWellPotential::laplacianV(const SingleVec& /* x */) {
    // @todo Complete the Laplacian?
    return 0.0;
}
