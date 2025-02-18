#include <array>
#include <fstream>
#include <cmath>
#include "mpi.h"
#include "bosonic_exchange_shuffle.h"
#include "simulation.h"

BosonicExchangeShuffle::BosonicExchangeShuffle(const Simulation& _sim) : BosonicExchange(_sim) {}

/**
 * @brief Shuffle temporary bosons order.
 */
void BosonicExchangeShuffle::assignIndirectionCoords() {
    if (sim.this_bead == 0) {
        sim.mars_gen->shuffle(indexes);
        MPI_Send(indexes.data(), nbosons, MPI_INT, nbeads - 1, 0, MPI_COMM_WORLD);
    } else {
        MPI_Status status;
        MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(indexes.data(), nbosons, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    for (int i = 0; i < nbosons; i++) {
        for (int axis = 0; axis < NDIM; ++axis) {
            indirection_x(i, axis) = x(indexes[i], axis);
            indirection_x_prev(i, axis) = x_prev(indexes[i], axis);
            indirection_x_next(i, axis) = x_next(indexes[i], axis);
        }
    }
}
