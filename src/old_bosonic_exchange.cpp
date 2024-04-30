#include <numeric>
#include <algorithm>
#include "mpi.h"

#include "old_bosonic_exchange.h"

OldBosonicExchange::OldBosonicExchange(int nbosons_, int np_, int bead_num_, double beta_, double spring_constant_,
                                       const dVec& x_, const dVec& x_prev_, const dVec& x_next_, bool pbc_,
                                       double size_) :
    BosonicExchangeBase(nbosons_, np_, bead_num_, beta_, spring_constant_, x_, x_prev_, x_next_, pbc_, size_),
    labels(nbosons_) {
    // Fill the labels array with numbers from 0 to nbosons-1
    std::iota(labels.begin(), labels.end(), 0);

    // For numerical stability
    e_shift = getMaxExteriorSpringEnergy();
}

/**
 * @brief Recalculates the longest exterior spring energy after each coordinate update.
 */
void OldBosonicExchange::prepare() {
    e_shift = getMaxExteriorSpringEnergy();
}

/**
 * Identifies the particle to which the neighboring P bead, connected to the first bead of ptcl_idx, belongs.
 * In the context of Cauchy's two-line notation, it retrieves the number from the upper line, placed above ptcl_idx.
 * Equivalently, it is the position of ptcl_idx in the bottom line, assuming zero-based indexing.
 * This works because the upper line can be viewed as the list of particles in the last imaginary time slice and
 * the bottom line as the list of their respective neighbors in the first imaginary time slice.
 *
 * @param ptcl_idx Index of the particle associated with the first bead.
 * @return Index of the particle associated with the neighboring P bead.
 */
int OldBosonicExchange::firstBeadNeighbor(int ptcl_idx) const {
    return std::distance(labels.begin(), std::ranges::find(labels, ptcl_idx));
}

/**
 * Identifies the particle to which the neighboring 1 bead, connected to the last bead of ptcl_idx, belongs.
 * In the context of Cauchy's two-line notation, it retrieves the number from the bottom line, placed below ptcl_idx.
 * Equivalently, it is the particle label in the bottom line at position ptcl_idx, assuming zero-based indexing.
 * This works because the upper line can be viewed as the list of particles in the last imaginary time slice and
 * the bottom line as the list of their respective neighbors in the first imaginary time slice.
 *
 * @param ptcl_idx Index of the particle associated with the last bead.
 * @return Index of the particle associated with the neighboring 1 bead.
 */
int OldBosonicExchange::lastBeadNeighbor(int ptcl_idx) const {
    return labels[ptcl_idx];
}

/**
 * Calculates the largest exterior spring energy that a permutation can yield in a given time-step. This energy is then used to shift
 * the spring energies in the Boltzmann weights to avoid numerical instabilities, which can be especially prominent
 * at high Trotter numbers.
 *
 * @return The largest exterior spring energy contribution due to a permutation.
 */
double OldBosonicExchange::getMaxExteriorSpringEnergy() {
    // Given the fact that the exterior spring energy is symmetric with respect to the first and the last time-slice,
    // we can calculate the exterior spring energy at the last time-slice and then communicate it to the first time-slice.
    if (bead_num == nbeads - 1) {
        double max_delta = 0.0;

        // Iterate over all permutations
        do {
            double delta_e_sigma = 0.0;

            // Iterate over all particles at time-slice P (within the current permutation)
            for (int l = 0; l < nbosons; ++l) {
                std::vector<double> sums(NDIM, 0.0);

                double diff_next[NDIM];

                // First bead (1) of some particle (depending on the permutation) minus last bead (P) of the l-th particle
                getBeadsSeparation(x, l, x_next, lastBeadNeighbor(l), diff_next);

                // Coordinate differences corresponding to exterior beads
                for (int axis = 0; axis < NDIM; ++axis) {
                    delta_e_sigma += 0.5 * spring_constant * diff_next[axis] * diff_next[axis];
                }
            }

            // Compare the current total exterior spring energy with the maximum total exterior spring energy found so far
            max_delta = std::max(max_delta, delta_e_sigma);
        } while (std::ranges::next_permutation(labels).found);

        // Communicate the maximum exterior spring energy to the first time-slice
        MPI_Send(&max_delta, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

        return max_delta;
    }

    if (bead_num == 0) {
        double max_delta_from_np;

        MPI_Recv(&max_delta_from_np, 1, MPI_DOUBLE, nbeads - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        return max_delta_from_np;
    }

    return 0.0;
}

/**
 * Calculates the contribution of the current time-slice to the effective bosonic spring potential.
 * For interior beads, this is the same as the ordinary spring energies in systems of distinguishable particles.
 * For exterior beads, the contribution is no longer harmonic, but is a sum of Boltzmann weights due to all permutations.
 *
 * @return Effective bosonic exchange potential.
 */
double OldBosonicExchange::effectivePotential() {
    // If the current time-slice is not the last one, then the spring energy is the 
    // same as the interior spring energy.
    if (bead_num != nbeads - 1)
        return interiorSpringEnergy();

    long permutation_counter = 0;
    double weights_sum = 0.0;

    // Iterate over all permutations and calculate the weights
    // associated with the exterior beads.
    do {
        permutation_counter++;
        double diff2 = 0.0;

        for (int ptcl_idx = 0; ptcl_idx < nbosons; ++ptcl_idx) {
            double diff[NDIM];

            getBeadsSeparation(x, ptcl_idx, x_next, lastBeadNeighbor(ptcl_idx), diff);

            // Coordinate differences corresponding to exterior beads
            for (int axis = 0; axis < NDIM; ++axis) {
                diff2 += diff[axis] * diff[axis];
            }
        }

        weights_sum += exp(-beta * 0.5 * spring_constant * diff2);
    } while (std::ranges::next_permutation(labels).found);

    /// @todo Think about how the i-Pi convention might affect the calculation
    return (-1.0 / beta) * log(weights_sum / permutation_counter);
}

/**
 * Evaluates the bosonic spring forces acting on the particles at the last imaginary time slice.
 *
 * @param[out] f Spring forces acting on the particles at time-slice P.
 */
void OldBosonicExchange::springForceLastBead(dVec& f) {
    const dVec x_first_bead = x_next;
    const dVec x_last_bead = x;

    dVec temp_force(nbosons);
    double denom_weight = 0.0;

    do {
        double weight = 0.0;

        for (int l = 0; l < nbosons; ++l) {
            std::vector<double> sums(NDIM, 0.0);

            double diff_next[NDIM];

            // First bead (1) of some particle (depending on the permutation) minus last bead (P) of the l-th particle
            getBeadsSeparation(x_last_bead, l, x_first_bead, lastBeadNeighbor(l), diff_next);

            // Coordinate differences corresponding to exterior beads
            for (int axis = 0; axis < NDIM; ++axis) {
                weight += diff_next[axis] * diff_next[axis];
            }

            double diff_prev[NDIM];

            // Previous bead (P-1) of the l-th particle minus last bead (P) of the l-th particle
            getBeadsSeparation(x_last_bead, l, x_prev, l, diff_prev);

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] += diff_prev[axis] + diff_next[axis];
                temp_force(l, axis) = sums[axis] * spring_constant;
            }
        }

        weight = exp(-beta * (0.5 * spring_constant * weight - e_shift));

        for (int l = 0; l < nbosons; ++l) {
            for (int axis = 0; axis < NDIM; ++axis) {
                f(l, axis) += weight * temp_force(l, axis);
            }
        }

        denom_weight += weight;
    } while (std::ranges::next_permutation(labels).found);

    // Normalize the forces by the sum of Boltzmann contributions due to all the permutations
    for (int l = 0; l < nbosons; ++l) {
        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) = f(l, axis) / denom_weight;
        }
    }
}

/**
 * Evaluates the bosonic spring forces acting on the particles at the first imaginary time slice.
 *
 * @param[out] f Spring forces acting on the particles at time-slice 1.
 */
void OldBosonicExchange::springForceFirstBead(dVec& f) {
    const dVec x_first_bead = x;
    const dVec x_last_bead = x_prev;

    dVec temp_force(nbosons);
    double denom_weight = 0.0;

    do {
        double weight = 0.0;

        for (int l = 0; l < nbosons; ++l) {
            std::vector<double> sums(NDIM, 0.0);

            double diff_prev[NDIM];

            // Last bead (P) of some particle (depending on the permutation) minus first bead (1) of the l-th particle
            getBeadsSeparation(x_first_bead, l, x_last_bead, firstBeadNeighbor(l), diff_prev);

            // Coordinate differences corresponding to exterior beads
            for (int axis = 0; axis < NDIM; ++axis) {
                weight += diff_prev[axis] * diff_prev[axis];
            }

            double diff_next[NDIM];

            // Next bead (2) of the l-th particle minus first bead (1) of the l-th particle
            getBeadsSeparation(x_first_bead, l, x_next, l, diff_next);

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] += diff_prev[axis] + diff_next[axis];
                temp_force(l, axis) = sums[axis] * spring_constant;
            }
        }

        weight = exp(-beta * (0.5 * spring_constant * weight - e_shift));

        for (int l = 0; l < nbosons; ++l) {
            for (int axis = 0; axis < NDIM; ++axis) {
                f(l, axis) += weight * temp_force(l, axis);
            }
        }

        denom_weight += weight;
    } while (std::ranges::next_permutation(labels).found);

    // Normalize the forces by the sum of Boltzmann contributions due to all the permutations
    for (int l = 0; l < nbosons; ++l) {
        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) = f(l, axis) / denom_weight;
        }
    }
}

/**
 * Calculates the contribution of the exterior beads to the primitive kinetic energy estimator.
 *
 * @return Weighted average of exterior spring energies over all permutations. 
 */
double OldBosonicExchange::primEstimator() {
    if (bead_num != 0)
        throw std::runtime_error("primEstimator() can only be called at the 1st time-slice");

    double numerator = 0.0;
    double denom_weight = 0.0;

    do {
        double weight = 0.0;

        for (int l = 0; l < nbosons; ++l) {
            double diff_prev[NDIM];

            // Last bead (P) of some particle (depending on the permutation) minus first bead (1) of the l-th particle
            getBeadsSeparation(x, l, x_prev, firstBeadNeighbor(l), diff_prev);

            // Coordinate differences corresponding to exterior beads
            for (int axis = 0; axis < NDIM; ++axis) {
                weight += diff_prev[axis] * diff_prev[axis];
            }
        }

        double delta_spring_energy = 0.5 * spring_constant * weight;
        weight = exp(-beta * (delta_spring_energy - e_shift));

        // Delta(E^sigma) is ready only at this point. We multiply Delta(E^sigma) by the weight
        // and add the result to the numerator (and also update the denominator by adding the 
        // weight of the current permutation).
        numerator += delta_spring_energy * weight;
        denom_weight += weight;
    } while (std::ranges::next_permutation(labels).found);

    double bosonic_spring_energy = numerator / denom_weight;

#if IPI_CONVENTION
    bosonic_spring_energy /= nbeads;
#endif

    return bosonic_spring_energy;
}
