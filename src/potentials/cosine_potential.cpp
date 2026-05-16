#include "potentials/cosine_potential.h"
#include <numbers>

CosinePotential::CosinePotential(double amplitude, double wavelength, double phase)
    : amplitude(amplitude), wavelength(wavelength), phase(phase) {
    k = 2 * std::numbers::pi / wavelength;
}

double CosinePotential::V(const SingleVec& x) {
    double potential = 0.0;
    for (int axis = 0; axis < NDIM; ++axis) {
        potential += std::cos(k * x[axis] + phase);
    }
    return potential * amplitude;
}

SingleVec CosinePotential::gradV(const SingleVec& x) {
    SingleVec grad{};
    const double prefactor = -amplitude * k;
    for (int axis = 0; axis < NDIM; ++axis) {
        grad[axis] = prefactor * std::sin(k * x[axis] + phase);
    }
    return grad;
}

double CosinePotential::laplacianV(const SingleVec& /* x */) {
    // @todo Complete the Laplacian?
    return 0.0;
}