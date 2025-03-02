#include <array>
#include <fstream>
#include <cmath>
#include "mpi.h"
#include <algorithm> 

#include "bosonic_exchange_PIS.h"
#include "simulation.h"

BosonicExchangePIS::BosonicExchangePIS(const Simulation& _sim) : BosonicExchange(_sim) {
    // current_measure = measure_shuffle(indexes);
}

/**
 * @brief Shuffle temporary bosons order.
 */
// void BosonicExchangePIS::assignIndirectionCoords() {
//     if (sim.this_bead == 0) {
//         for (int i = 0; i < 1000; i++) {
//             std::vector<int>& suggested_indexes = indexes;
//             sim.mars_gen->shuffle(suggested_indexes);
//             double new_measure = measure_shuffle(suggested_indexes);
//             // double random_uniform = sim.mars_gen->uniform();
//             // if (random_uniform < current_measure / (current_measure + new_measure)) {
//             //     indexes = suggested_indexes;
//             //     current_measure = new_measure;
//             // }
//             if (new_measure < current_measure) {
//                 indexes = suggested_indexes;
//                 current_measure = new_measure;
//             }
//         }
//         printStatus(std::format("{:.3}", current_measure), 0);
//         MPI_Send(indexes.data(), nbosons, MPI_INT, nbeads - 1, 0, MPI_COMM_WORLD);
//     } else {
//         MPI_Status status;
//         MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
//         MPI_Recv(indexes.data(), nbosons, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//     }
    
//     for (int i = 0; i < nbosons; i++) {
//         for (int axis = 0; axis < NDIM; ++axis) {
//             indirection_x(i, axis) = x(indexes[i], axis);
//             indirection_x_prev(i, axis) = x_prev(indexes[i], axis);
//             indirection_x_next(i, axis) = x_next(indexes[i], axis);
//         }
//     }
// }

// double BosonicExchangePIS::measure_shuffle(std::vector<int>& suggested_indexes) const {
//     double new_measure = getBeadsSeparationSquared(x, suggested_indexes[0], x_prev, suggested_indexes[nbosons - 1]);
//     for (int i = 1; i < nbosons; i++) {
//         new_measure += getBeadsSeparationSquared(x, suggested_indexes[i], x_prev, suggested_indexes[i - 1]);
//     }
//     return new_measure;
// }

void BosonicExchangePIS::assignIndirectionCoords() {
    if (sim.this_bead == 0) {
        std::vector<int> suggested_indexes(indexes.begin() + 1, indexes.end());
        int current_index = indexes[0];
        for (int i = 0; i < nbosons-1; i++) {
            current_index = get_next_index(current_index, suggested_indexes);
            indexes[i+1] = current_index;
            suggested_indexes.erase(std::remove(suggested_indexes.begin(), suggested_indexes.end(), current_index), suggested_indexes.end());
            // printStatus(std::format("{}", current_index), 0);
        }
        // printStatus("\n \n \n \n", 0);
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

int BosonicExchangePIS::get_next_index(int current_index, std::vector<int>& suggested_indexes) const {
    double current_measure = getBeadsSeparationSquared(x, current_index, x_next, suggested_indexes[0]);
    // printStatus(std::format("current, {:.3}", current_measure), 0);
    int next_index = suggested_indexes[0];
    for (int i = 1; i < suggested_indexes.size(); i++) {
        double new_measure = getBeadsSeparationSquared(x, current_index, x_next, suggested_indexes[i]);
        // printStatus(std::format("{:.3}", new_measure), 0);
        if (new_measure < current_measure) {
            // printStatus(std::format("accepted, {}", suggested_indexes[i]), 0);
            current_measure = new_measure;
            next_index = suggested_indexes[i];
        }
    }
    return next_index;
}