#include "old_bosonic_exchange.h"
#include "mpi.h"

#include <numeric>

OldBosonicExchange::OldBosonicExchange(int nbosons, int np, int bead_num, double _beta, double spring_constant,
    const dVec x, const dVec x_prev, const dVec x_next, bool pbc, double L) : 
    BosonicExchangeBase(nbosons, np, bead_num, _beta, spring_constant, x, x_prev, x_next, pbc, L),
    labels(nbosons)
{
    // Fill the labels array with numbers from 0 to nbosons-1
    std::iota(labels.begin(), labels.end(), 0);

    // For numerical stability
    e_longest = get_elongest();
}

/* ---------------------------------------------------------------------- */

OldBosonicExchange::~OldBosonicExchange() 
{
}

/* ---------------------------------------------------------------------- */

void OldBosonicExchange::updateCoordinates(const dVec new_x, const dVec new_x_prev, const dVec new_x_next)
{
    BosonicExchangeBase::updateCoordinates(new_x, new_x_prev, new_x_next);

    // Recalculate e_longest after each coordinate update
    e_longest = get_elongest();
}

/* ---------------------------------------------------------------------- */

// Tells which particle the neighboring P bead (which is connected to 1st bead of ptcl_idx) belongs to
int OldBosonicExchange::neighbor_of_first(int ptcl_idx)
{
    auto it = std::find(labels.begin(), labels.end(), ptcl_idx);

    return std::distance(labels.begin(), it);
}

/* ---------------------------------------------------------------------- */

// Tells which particle the neighboring 1 bead (which is connected to P bead of ptcl_idx) belongs to
int OldBosonicExchange::neighbor_of_last(int ptcl_idx)
{
    return labels[ptcl_idx];
}

/* ---------------------------------------------------------------------- */

double OldBosonicExchange::get_elongest()
{
    const dVec x_first_bead = x_next;
    const dVec x_last_bead = x;

    if (bead_num == np - 1) {
        double max_delta = 0.0;

        std::vector<int> local_labels(nbosons);
        std::iota(local_labels.begin(), local_labels.end(), 0);

        // Iterate over all permutations
        do {
            double delta_e_sigma = 0.0;

            // Iterate over all particles
            for (int l = 0; l < nbosons; ++l) {
                std::vector<double> sums(NDIM, 0.0);

                double diff_next[NDIM];

                // First bead (1) of some particle (depending on the permutation) minus last bead (P) of the l-th particle
                diff_two_beads(x_last_bead, l, x_first_bead, neighbor_of_last(l), diff_next);

                // Coordinate differences corresponding to exterior beads
                for (int axis = 0; axis < NDIM; ++axis) {
                    delta_e_sigma += 0.5 * spring_constant * diff_next[axis] * diff_next[axis];
                }
            }

            max_delta = std::max(max_delta, delta_e_sigma);
        } while (std::next_permutation(local_labels.begin(), local_labels.end()));

        MPI_Send(&max_delta, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

        return max_delta;

    }
    else if (bead_num == 0) {
        double max_delta_from_np;

        MPI_Recv(&max_delta_from_np, 1, MPI_DOUBLE, np - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        return max_delta_from_np;
    }

    return 0.0;
}

/* ---------------------------------------------------------------------- */

void OldBosonicExchange::spring_force_last_bead(dVec& f)
{
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
            diff_two_beads(x_last_bead, l, x_first_bead, neighbor_of_last(l), diff_next);

            // Coordinate differences corresponding to exterior beads
            for (int axis = 0; axis < NDIM; ++axis) {
                weight += diff_next[axis] * diff_next[axis];
            }

            double diff_prev[NDIM];

            // Previous bead (P-1) of the l-th particle minus last bead (P) of the l-th particle
            diff_two_beads(x_last_bead, l, x_prev, l, diff_prev);

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] += diff_prev[axis] + diff_next[axis];
                temp_force(l, axis) = sums[axis] * spring_constant;
            }
        }

        weight = exp(-beta * (0.5 * spring_constant * weight - e_longest));

        for (int l = 0; l < nbosons; ++l) {
            for (int axis = 0; axis < NDIM; ++axis) {
                f(l, axis) += weight * temp_force(l, axis);
            }
        }

        denom_weight += weight;

    } while (std::next_permutation(labels.begin(), labels.end()));

    // Normalize the forces by the sum of Boltzmann contributions due to all the permutations
    for (int l = 0; l < nbosons; ++l) {
        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) = f(l, axis) / denom_weight;
        }
    }

    // Return labels to the original state (identity permutation)
    std::sort(labels.begin(), labels.end());
}

/* ---------------------------------------------------------------------- */

void OldBosonicExchange::spring_force_first_bead(dVec& f)
{
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
            diff_two_beads(x_first_bead, l, x_last_bead, neighbor_of_first(l), diff_prev);

            // Coordinate differences corresponding to exterior beads
            for (int axis = 0; axis < NDIM; ++axis) {
                weight += diff_prev[axis] * diff_prev[axis];
            }

            double diff_next[NDIM];

            // Next bead (2) of the l-th particle minus first bead (1) of the l-th particle
            diff_two_beads(x_first_bead, l, x_next, l, diff_next);

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] += diff_prev[axis] + diff_next[axis];
                temp_force(l, axis) = sums[axis] * spring_constant;
            }
        }

        weight = exp(-beta * (0.5 * spring_constant * weight - e_longest));

        for (int l = 0; l < nbosons; ++l) {
            for (int axis = 0; axis < NDIM; ++axis) {
                f(l, axis) += weight * temp_force(l, axis);
            }
        }

        denom_weight += weight;

    } while (std::next_permutation(labels.begin(), labels.end()));

    // Normalize the forces by the sum of Boltzmann contributions due to all the permutations
    for (int l = 0; l < nbosons; ++l) {
        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) = f(l, axis) / denom_weight;
        }
    }

    // Return labels to the original state (identity permutation)
    std::sort(labels.begin(), labels.end());
}

/* ---------------------------------------------------------------------- */

// This method should be called only at timeslice=0
double OldBosonicExchange::prim_estimator() {
    const dVec x_first_bead = x;
    const dVec x_last_bead = x_prev;

    double numerator = 0.0;
    double denom_weight = 0.0;

    do {
        double weight = 0.0;

        for (int l = 0; l < nbosons; ++l) {
            std::vector<double> sums(NDIM, 0.0);

            double diff_prev[NDIM];

            // Last bead (P) of some particle (depending on the permutation) minus first bead (1) of the l-th particle
            diff_two_beads(x_first_bead, l, x_last_bead, neighbor_of_first(l), diff_prev);

            // Coordinate differences corresponding to exterior beads
            for (int axis = 0; axis < NDIM; ++axis) {
                weight += diff_prev[axis] * diff_prev[axis];
            }
        }

        double delta_spring_energy = 0.5 * spring_constant * weight - e_longest;
        weight = exp(-beta * delta_spring_energy);

        // Delta(E^sigma) is ready only at this point. We multiply Delta(E^sigma) by the weight
        // and add the result to the numerator (and also update the denominator by adding the 
        // weight of the current permutation).
        numerator += delta_spring_energy * weight;
        denom_weight += weight;

    } while (std::next_permutation(labels.begin(), labels.end()));

    // Return labels to the original state (identity permutation)
    std::sort(labels.begin(), labels.end());

    double bosonic_spring_energy = numerator / denom_weight;

#if IPI_CONVENTION
    bosonic_spring_energy /= np;
#endif

    return bosonic_spring_energy;
}