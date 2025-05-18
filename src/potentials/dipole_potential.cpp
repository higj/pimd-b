#include "potentials/dipole_potential.h"

DipolePotential::DipolePotential(double strength) : strength(strength) {}

double DipolePotential::V(const dVec& x) {
    double potential = 0.0;

    for (int ptcl_idx = 0; ptcl_idx < x.len(); ++ptcl_idx) {
        /*
        double norm_squared = 0.0;

        for (int axis = 0; axis < NDIM; ++axis) {
            norm_squared += x(ptcl_idx, axis) * x(ptcl_idx, axis);
        }
        potential += strength / (norm_squared * sqrt(norm_squared));
        */
        const double norm = x.norm(ptcl_idx);
        potential += strength / (norm * norm * norm);
    }

    return potential;
}

dVec DipolePotential::gradV(const dVec& x) {
    dVec tempr(x);

    for (int ptcl_idx = 0; ptcl_idx < x.len(); ++ptcl_idx) {
        /*
        double norm_squared = 0.0;
        for (int axis = 0; axis < NDIM; ++axis) {
            norm_squared += x(ptcl_idx, axis) * x(ptcl_idx, axis);
        }

        const double norm = sqrt(norm_squared);
        */
        const double norm = x.norm(ptcl_idx);
        const double prefactor = -3.0 * strength / (norm * norm * norm * norm * norm);

        for (int axis = 0; axis < NDIM; ++axis) {
            tempr(ptcl_idx, axis) = prefactor * tempr(ptcl_idx, axis);
        }
    }

    return tempr;
}

double DipolePotential::laplacianV(const dVec& x) {
    // @todo Complete the Laplacian?
    return 0.0;
}