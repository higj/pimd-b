#include <array>
#include <fstream>
#include <cmath>

#include "bosonic_exchange_carrousel.h"
#include "simulation.h"

BosonicExchangeCarrousel::BosonicExchangeCarrousel(const Simulation& _sim) : BosonicExchange(_sim), 
    temp_x(_sim.natoms), temp_x_prev(_sim.natoms), temp_x_next(_sim.natoms) {}

/**
 * @brief Re-evaluate the bosonic energies and connection probabilities.
 * Typically used after coordinate updates.
 */
void BosonicExchangeCarrousel::prepare() {
    assignTempCoords();
    evaluateBosonicEnergies();
}

/**
 * @brief Set temporary bosons order accurding to current timestep.
 */
void BosonicExchangeCarrousel::assignTempCoords() {
    for (int i = 0; i < nbosons; i++) {
        int index = ((i + sim.getStep()) % nbosons);
        for (int axis = 0; axis < NDIM; ++axis) {
            temp_x(i, axis) = x(index, axis);
            temp_x_prev(i, axis) = x_prev(index, axis);
            temp_x_next(i, axis) = x_next(index, axis);
        }
    }
}

/**
 * Calculates the bosonic force exerted on the beads
 * at the first and the last time-slices.
 *
 * @param[out] f Vector to store the forces.
 */
void BosonicExchangeCarrousel::exteriorSpringForce(dVec& f) {
    dVec temp_f(nbosons);
    std::vector<int> indexes(nbosons);
    for (int i = 0; i < nbosons; i++) {
        indexes[i] = ((i + sim.getStep()) % nbosons);
        for (int axis = 0; axis < NDIM; ++axis) {
            temp_f(i, axis) = f(indexes[i], axis);
        }
    }

    if (sim.this_bead == 0) {
        springForceFirstBead(temp_f, temp_x, temp_x_prev, temp_x_next);
    } else {
        springForceLastBead(temp_f, temp_x, temp_x_prev, temp_x_next);
    }

    for (int i = 0; i < nbosons; i++) {
        for (int axis = 0; axis < NDIM; ++axis) {
            f(indexes[i], axis) = temp_f(i, axis);
        }
    }
}

/**
 * Determines the coordinates of the first and the last bead based on the current time-slice.
 *
 * @param[out] x_first_bead The coordinates of the particles in the first time-slice.
 * @param[out] x_last_bead The coordinates of the particles in the last time-slice.
 */
void BosonicExchangeCarrousel::assignFirstLast(dVec& x_first_bead, dVec& x_last_bead) const {
    if (sim.this_bead == 0) {
        x_first_bead = temp_x;
        x_last_bead = temp_x_prev;
    } else {
        x_first_bead = temp_x_next;
        x_last_bead = temp_x;
    }
}
