#include "potentials/anharmonic_potential.h"

AnharmonicPotential::AnharmonicPotential(double mass, double omega, double cubic_const, double quart_const)
    : mass(mass), omega(omega), cubic_const(cubic_const), quart_const(quart_const) {
    k = mass * omega * omega;
}

double AnharmonicPotential::V(const Vec& x) {
    double potential = 0.0;
    for (int axis = 0; axis < NDIM; ++axis) {
        double xi = x[axis];
        double xi2 = xi * xi;
        double xi3 = xi2 * xi;
        double xi4 = xi3 * xi;
        potential += 0.5 * k * xi2 + cubic_const * xi3 + quart_const * xi4;
    }
    return potential;
}

Vec AnharmonicPotential::gradV(const Vec& x) {
    Vec grad{};
    for (int axis = 0; axis < NDIM; ++axis) {
        double xi = x[axis];
        double xi2 = xi * xi;
        double xi3 = xi2 * xi;
        grad[axis] = k * xi + 3.0 * cubic_const * xi2 + 4.0 * quart_const * xi3;
    }
    return grad;
}

double AnharmonicPotential::laplacianV(const Vec& x) {
    double lap = 0.0;
    for (int axis = 0; axis < NDIM; ++axis) {
        double xi = x[axis];
        double xi2 = xi * xi;
        lap += k + 6.0 * cubic_const * xi + 12.0 * quart_const * xi2;
    }
    return lap;
}
