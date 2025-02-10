#include <array>
#include <fstream>
#include <cmath>

#include "bosonic_exchange_carrousel.h"
#include "simulation.h"

BosonicExchangeCarrousel::BosonicExchangeCarrousel(const Simulation& _sim) : BosonicExchange(_sim) {}

/**
 * @brief Set temporary bosons order accurding to current timestep.
 */
void BosonicExchangeCarrousel::assignIndirectionCoords() {
    for (int i = 0; i < nbosons; i++) {
        indexes[i] = ((i + sim.getStep()) % nbosons);
        for (int axis = 0; axis < NDIM; ++axis) {
            indirection_x(i, axis) = x(indexes[i], axis);
            indirection_x_prev(i, axis) = x_prev(indexes[i], axis);
            indirection_x_next(i, axis) = x_next(indexes[i], axis);
        }
    }
}
