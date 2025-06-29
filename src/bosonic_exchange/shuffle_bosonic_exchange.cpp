#include "common.h"
#include "bosonic_exchange/shuffle_bosonic_exchange.h"
#include "mpi.h"
BosonicExchangeShuffle::BosonicExchangeShuffle(Params& param_obj, const dVec& coord, const dVec& prev_coord, const dVec& next_coord, const int this_bead) :
    BosonicExchange(param_obj, coord, prev_coord, next_coord, this_bead), step(0), indexes(nbosons) {
    getVariant(param_obj.sim["shuffle_timer"], timer);	
    for (int i = 0; i < nbosons; ++i) {
        indexes[i] = i;
    }
    unsigned int params_seed;
    getVariant(param_obj.sim["seed"], params_seed);
    mars_gen = std::make_unique<RanMars>(params_seed);

}

void BosonicExchangeShuffle::assignIndirectionCoords() {
    if (step == timer) {
        if (this_bead == 0) {
             mars_gen->shuffle(indexes);
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
	step = 0;
    } else {
	step += 1;
    }
}
