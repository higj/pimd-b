#include "potential.h"

Potential::Potential() : tailV(0.0) {

}

Potential::~Potential() {
}

// Isotropic harmonic potential
HarmonicPotential::HarmonicPotential(double mass, double omega) : mass(mass), omega(omega) {
    k = mass * omega * omega;
}

HarmonicPotential::~HarmonicPotential() {
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
    dVec temp = x;
    return k * temp;
}

double HarmonicPotential::laplacianV(const dVec& x) {
    // TODO: Complete the Laplacian?
    return 0.0;
}

// Double-well potential
DoubleWellPotential::DoubleWellPotential(double mass, double strength, double loc) : mass(mass), strength(strength), loc(loc) {
}

DoubleWellPotential::~DoubleWellPotential() {
}

double DoubleWellPotential::V(const dVec& x) {
    double potential = 0.0;
    double loc2 = loc * loc;

    for (int ptcl_idx = 0; ptcl_idx < x.len(); ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            potential += (x(ptcl_idx, axis) * x(ptcl_idx, axis) - loc2) * (x(ptcl_idx, axis) * x(ptcl_idx, axis) - loc2);
        }
    }

    potential *= mass * strength;

    return potential;
}

dVec DoubleWellPotential::gradV(const dVec& x) {
    double loc2 = loc * loc;
    dVec tempr;
    tempr = x;

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
    // TODO: Complete the Laplacian?
    return 0.0;
}


// Dipole potential
DipolePotential::DipolePotential(double strength) : strength(strength) {
}

DipolePotential::~DipolePotential() {
}

double DipolePotential::V(const dVec& x) {
    double potential = 0.0;

    for (int ptcl_idx = 0; ptcl_idx < x.len(); ++ptcl_idx) {
        double norm_squared = 0.0;

        for (int axis = 0; axis < NDIM; ++axis) {
            norm_squared += x(ptcl_idx, axis) * x(ptcl_idx, axis);
        }

        potential += strength / (norm_squared * sqrt(norm_squared));
    }

    return potential;
}

dVec DipolePotential::gradV(const dVec& x) {
    dVec tempr;
    tempr = x;

    for (int ptcl_idx = 0; ptcl_idx < x.len(); ++ptcl_idx) {
        double norm_squared = 0.0;
        for (int axis = 0; axis < NDIM; ++axis) {
            norm_squared += x(ptcl_idx, axis) * x(ptcl_idx, axis);
        }

        double norm = sqrt(norm_squared);
        double prefactor = -3.0 * strength / (norm * norm * norm * norm * norm);

        for (int axis = 0; axis < NDIM; ++axis) {
            tempr(ptcl_idx, axis) = prefactor * tempr(ptcl_idx, axis);
        }
    }

    return tempr;
}

double DipolePotential::laplacianV(const dVec& x) {
    // TODO: Complete the Laplacian?
    return 0.0;
}