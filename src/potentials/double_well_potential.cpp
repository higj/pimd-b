#include "potentials/double_well_potential.h"

DoubleWellPotential::DoubleWellPotential(double mass, double strength, double loc)
    : mass(mass), strength(strength), loc(loc) {}

double DoubleWellPotential::V(const Vec& x) {
    double potential = 0.0;
    const double loc2 = loc * loc;

    for (int axis = 0; axis < NDIM; ++axis) {
        potential += (x[axis] * x[axis] - loc2) * (x[axis] * x[axis] - loc2);
    }

    return potential * mass * strength;
}

Vec DoubleWellPotential::gradV(const Vec& x) {
    Vec grad{};
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

double DoubleWellPotential::laplacianV(const Vec& /* x */) {
    // @todo Complete the Laplacian?
    return 0.0;
}
