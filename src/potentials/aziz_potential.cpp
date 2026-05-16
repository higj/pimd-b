#include "potentials/aziz_potential.h"

AzizPotential::AzizPotential() {
    // Aziz potential (HFDHE2) parameters, based on [J. Chem. Phys. 70, 4330-4342 (1979)].
    // The dimensional quantities (rm and epsilon) are in atomic units.
    rm = 5.60738;            // 5.60738 Bohr = 2.9673 Angstrom
    A = 0.5448504e6;
    epsilon = 3.42016E-5;    // 3.42016E-5 Hartrees = kB * (10.8 Kelvin)
    alpha = 13.353384;
    D = 1.241314;
    C6 = 1.3732412;
    C8 = 0.4253785;
    C10 = 0.1781;

    /// @todo Compute the tail correction here
}

double AzizPotential::V(const SingleVec& x) {
    // Scale the distance relative to the Aziz equilibrium distance
    const double r_scaled = norm(x) / rm;
    const double Urep = A * exp(-alpha * r_scaled);

    if (r_scaled > EPS && r_scaled < 0.01) {
        // If the distance is small, the 6-8-10 terms are negligible
        // and the repulsion is dominated by "Urep"
        return epsilon * Urep;
    } else {
        const double ix2 = 1.0 / (r_scaled * r_scaled);
        const double ix6 = ix2 * ix2 * ix2;
        const double ix8 = ix6 * ix2;
        const double ix10 = ix8 * ix2;
        const double Uatt = -(C6 * ix6 + C8 * ix8 + C10 * ix10) * F(r_scaled);
        return epsilon * (Urep + Uatt);
    }
}

SingleVec AzizPotential::gradV(const SingleVec& x) {
    SingleVec grad{};
    const double n = norm(x);
    const double r_scaled = n / rm;
    const double T1 = -A * alpha * exp(-alpha * r_scaled);

    double grad_result;

    if (r_scaled > EPS && r_scaled < 0.01) {
        grad_result = T1 * (epsilon / rm);
    } else {
        const double ix = 1.0 / r_scaled;
        const double ix2 = ix * ix;
        const double ix6 = ix2 * ix2 * ix2;
        const double ix7 = ix6 * ix;
        const double ix8 = ix6 * ix2;
        const double ix9 = ix8 * ix;
        const double ix10 = ix8 * ix2;
        const double ix11 = ix10 * ix;

        const double T2 = (6.0 * C6 * ix7 + 8.0 * C8 * ix9 + 10.0 * C10 * ix11) * F(r_scaled);
        const double T3 = -(C6 * ix6 + C8 * ix8 + C10 * ix10) * dF(r_scaled);

        grad_result = (epsilon / rm) * (T1 + T2 + T3);
    }

    grad_result /= n;

    for (int axis = 0; axis < NDIM; ++axis) {
        grad[axis] = grad_result * x[axis];
    }

    return grad;
}

double AzizPotential::laplacianV(const SingleVec& x) {
    const double n = norm(x);
    const double r_scaled = n / rm;

    const double T1 = A * alpha * alpha * exp(-alpha * r_scaled);

    if (r_scaled > EPS && r_scaled < 0.01) {
        return (epsilon / rm) * T1;
    } else {
        const double ix = 1.0 / r_scaled;
        const double ix2 = ix * ix;
        const double ix6 = ix2 * ix2 * ix2;
        const double ix7 = ix6 * ix;
        const double ix8 = ix6 * ix2;
        const double ix9 = ix8 * ix;
        const double ix10 = ix8 * ix2;
        const double ix11 = ix10 * ix;
        const double ix12 = ix11 * ix;
        const double T2 = -(42.0 * C6 * ix8 + 72.0 * C8 * ix10 + 110.0 * C10 * ix12) * F(r_scaled);
        const double T3 = 2.0 * (6.0 * C6 * ix7 + 8.0 * C8 * ix9 + 10.0 * C10 * ix11) * dF(r_scaled);
        const double T4 = -(C6 * ix6 + C8 * ix8 + C10 * ix10) * d2F(r_scaled);
        return (epsilon / (rm * rm)) * (T1 + T2 + T3 + T4);
    }
}