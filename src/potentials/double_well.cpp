#include "potentials/double_well.h"

DoubleWellPotential::DoubleWellPotential(double mass, double strength, double loc)
    : mass(mass), strength(strength), loc(loc) {}

double DoubleWellPotential::V(const dVec& x) {
    double potential = 0.0;
    const double loc2 = loc * loc;

    for (int ptcl_idx = 0; ptcl_idx < x.len(); ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            potential += (x(ptcl_idx, axis) * x(ptcl_idx, axis) - loc2) *
                (x(ptcl_idx, axis) * x(ptcl_idx, axis) - loc2);
        }
    }

    potential *= mass * strength;

    return potential;
}

dVec DoubleWellPotential::gradV(const dVec& x) {
    dVec tempr(x);
    const double loc2 = loc * loc;

    for (int ptcl_idx = 0; ptcl_idx < x.len(); ++ptcl_idx) {
        double prefactor = 0;
        for (int axis = 0; axis < NDIM; ++axis) {
            prefactor += x(ptcl_idx, axis) * x(ptcl_idx, axis);
        }

        prefactor = 4 * mass * strength * (prefactor - loc2);

        for (int axis = 0; axis < NDIM; ++axis) {
            tempr(ptcl_idx, axis) = prefactor * tempr(ptcl_idx, axis);
        }
    }

    return tempr;
}

double DoubleWellPotential::laplacianV(const dVec& x) {
    // @todo Complete the Laplacian?
    return 0.0;
}
