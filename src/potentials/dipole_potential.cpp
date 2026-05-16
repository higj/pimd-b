#include "potentials/dipole_potential.h"

DipolePotential::DipolePotential(double strength) : strength(strength) {}

double DipolePotential::V(const Vec& x) {
    const double n = norm(x);
    return strength / (n * n * n);
}

Vec DipolePotential::gradV(const Vec& x) {
    Vec grad{};
    const double n = norm(x);
    const double prefactor = -3.0 * strength / (n * n * n * n * n);

    for (int axis = 0; axis < NDIM; ++axis) {
        grad[axis] = prefactor * x[axis];
    }
    return grad;
}

double DipolePotential::laplacianV(const Vec& /* x */) {
    // @todo Complete the Laplacian?
    return 0.0;
}