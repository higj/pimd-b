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
    shuffle();
}

void BosonicExchangeShuffle::assignIndirectionCoords() {
    if (step == timer) {
	shuffle();
	step = 0;
    } else {
	step += 1;
    }
    for (int i = 0; i < nbosons; i++) {
        for (int axis = 0; axis < NDIM; ++axis) {
            indirection_x(i, axis) = x(indexes[i], axis);
            indirection_x_prev(i, axis) = x_prev(indexes[i], axis);
            indirection_x_next(i, axis) = x_next(indexes[i], axis);
        }
    }
}

void BosonicExchangeShuffle::shuffle() {
    if (this_bead == 0) {
         mars_gen->shuffle(indexes);
         MPI_Send(indexes.data(), nbosons, MPI_INT, nbeads - 1, 0, MPI_COMM_WORLD);
    } else {
       MPI_Status status;
       MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
       MPI_Recv(indexes.data(), nbosons, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void BosonicExchangeShuffle::springForceLastBead(dVec& f) {
    for (int l = 0; l < nbosons; l++) {
        std::array<double, NDIM> sums = {};

        for (int next_l = 0; next_l <= l + 1 && next_l < nbosons; next_l++) {
            double diff_next[NDIM];

            getBeadsSeparation(indirection_x, l, indirection_x_next, next_l, diff_next);

            double prob = connection_probabilities[nbosons * l + next_l];

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] += prob * diff_next[axis];
            }
        }

        double diff_prev[NDIM];
        getBeadsSeparation(indirection_x, l, indirection_x_prev, l, diff_prev);

        for (int axis = 0; axis < NDIM; ++axis) {
            sums[axis] += diff_prev[axis];
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            f(indexes[l], axis) = sums[axis] * spring_constant;
        }
    }
}

void BosonicExchangeShuffle::springForceFirstBead(dVec& f) {
    for (int l = 0; l < nbosons; l++) {
        std::array<double, NDIM> sums = {};

        for (int prev_l = std::max(0, l - 1); prev_l < nbosons; prev_l++) {
            double diff_prev[NDIM];

            getBeadsSeparation(indirection_x, l, indirection_x_prev, prev_l, diff_prev);

            double prob = connection_probabilities[nbosons * prev_l + l];

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] += prob * diff_prev[axis];
            }
        }

        double diff_next[NDIM];
        getBeadsSeparation(indirection_x, l, indirection_x_next, l, diff_next);

        for (int axis = 0; axis < NDIM; ++axis) {
            sums[axis] += diff_next[axis];
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            f(indexes[l], axis) = sums[axis] * spring_constant;
        }
    }
}
