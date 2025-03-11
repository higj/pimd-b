#include "potentials/cosine.h"
#include <numbers>

CosinePotential::CosinePotential(double amplitude, double wavelength, double phase)
    : amplitude(amplitude), wavelength(wavelength), phase(phase) {
    k = 2 * std::numbers::pi / wavelength;
}

double CosinePotential::V(const dVec& x) {
    double potential = 0.0;

    for (int ptcl_idx = 0; ptcl_idx < x.len(); ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            potential += std::cos(k * x(ptcl_idx, axis) + phase);
        }
    }

    potential *= amplitude;

    return potential;
}

dVec CosinePotential::gradV(const dVec& x) {
    dVec tempr(x);

    const double prefactor = -amplitude * k;

    for (int ptcl_idx = 0; ptcl_idx < x.len(); ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            tempr(ptcl_idx, axis) = prefactor * std::sin(k * tempr(ptcl_idx, axis) + phase);
        }
    }

    return tempr;
}

double CosinePotential::laplacianV(const dVec& x) {
    // @todo Complete the Laplacian?
    return 0.0;
}