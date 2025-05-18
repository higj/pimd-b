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

double AzizPotential::V(const dVec& x) {
    double potential = 0.0;

    // Iterate through the N provided distance vectors
    for (int ptcl_idx = 0; ptcl_idx < x.len(); ++ptcl_idx) {
        // Calculate the scalar distance of a given pair of particles
        /*
        double norm_squared = 0.0;

        for (int axis = 0; axis < NDIM; ++axis) {
            norm_squared += x(ptcl_idx, axis) * x(ptcl_idx, axis);
        }

        const double r_scaled = sqrt(norm_squared) / rm;
        */

        // Scale the distance relative to the Aziz equilibrium distance
        const double r_scaled = x.norm(ptcl_idx) / rm;

        const double Urep = A * exp(-alpha * r_scaled);

        if (r_scaled > EPS && r_scaled < 0.01) {
            // If the distance is small, the 6-8-10 terms are negligible
            // and the repulsion is dominated by "Urep"
            potential += epsilon * Urep;
        } else {
            const double ix2 = 1.0 / (r_scaled * r_scaled);
            const double ix6 = ix2 * ix2 * ix2;
            const double ix8 = ix6 * ix2;
            const double ix10 = ix8 * ix2;
            const double Uatt = -(C6 * ix6 + C8 * ix8 + C10 * ix10) * F(r_scaled);
            potential += epsilon * (Urep + Uatt);
        }
    }

    return potential;
}

dVec AzizPotential::gradV(const dVec& x) {
    dVec tempr(x);

    // Iterate through the N provided distance vectors
    for (int pair_idx = 0; pair_idx < x.len(); ++pair_idx) {
        // Calculate the scalar distance of a given pair of particles
        /*
        double norm_squared = 0.0;

        for (int axis = 0; axis < NDIM; ++axis) {
            norm_squared += x(pair_idx, axis) * x(pair_idx, axis);
        }

        const double norm = sqrt(norm_squared);
        */
        const double norm = x.norm(pair_idx);
        const double r_scaled = norm / rm;
        const double T1 = -A * alpha * exp(-alpha * r_scaled);

        double grad_result;

        // Do not allow self-interactions
        if (r_scaled > EPS && r_scaled < 0.01)
            grad_result = T1 * (epsilon / rm);
        else {
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

        // The gradient multiplies the distance vector.
        // We need only the direction, so we normalize the distance vector.
        grad_result /= norm;

        for (int axis = 0; axis < NDIM; ++axis) {
            tempr(pair_idx, axis) = grad_result * tempr(pair_idx, axis);
        }
    }

    return tempr;
}

double AzizPotential::laplacianV(const dVec& x) {
    double laplacian = 0.0;

    // Iterate through the N provided distance vectors
    for (int pair_idx = 0; pair_idx < x.len(); ++pair_idx) {
        // Calculate the scalar distance of a given pair of particles
        /*
        double norm_squared = 0.0;

        for (int axis = 0; axis < NDIM; ++axis) {
            norm_squared += x(pair_idx, axis) * x(pair_idx, axis);
        }

        const double norm = sqrt(norm_squared);
        */
        const double norm = x.norm(pair_idx);
        const double r_scaled = norm / rm;

        const double T1 = A * alpha * alpha * exp(-alpha * r_scaled);

        // Hard core limit
        if (r_scaled > EPS && r_scaled < 0.01) {
            laplacian += (epsilon / rm) * T1;
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
            laplacian += (epsilon / (rm * rm)) * (T1 + T2 + T3 + T4);
        }
    }

    return laplacian;
}