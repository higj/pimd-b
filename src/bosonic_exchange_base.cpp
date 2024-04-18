#include "bosonic_exchange_base.h"

BosonicExchangeBase::BosonicExchangeBase(int nbosons_, int np_, int bead_num_, double beta_, double spring_constant_,
                                         const dVec& x_, const dVec& x_prev_, const dVec& x_next_, bool pbc_,
                                         double size_) :
    nbosons(nbosons_),
    nbeads(np_),
    bead_num(bead_num_),
    spring_constant(spring_constant_),
    beta(beta_),
    x(x_),
    x_prev(x_prev_),
    x_next(x_next_),
    pbc(pbc_),
    size(size_) {
#if IPI_CONVENTION
    beta = beta_ / nbeads;
#endif
}

/* ---------------------------------------------------------------------- */

void BosonicExchangeBase::getBeadsSeparation(const dVec x1, int l1, const dVec x2, int l2, double diff[NDIM]) const {
    l1 = l1 % nbosons;
    l2 = l2 % nbosons;

    for (int axis = 0; axis < NDIM; ++axis) {
        double dx = x2(l2, axis) - x1(l1, axis);
#if MINIM
        if (pbc)
            applyMinimumImage(dx, size);
#endif
        diff[axis] = dx;
    }
}

/**
 * Calculates the distance squared between two beads.
 * 
 * @param x1 Coordinates of the first bead
 * @param l1 Particle index of the first bead
 * @param x2 Coordinates of the second bead
 * @param l2 Particle index of the second bead
 * @return Distance squared between the two beads
 */
double BosonicExchangeBase::getBeadsSeparationSquared(const dVec x1, int l1, const dVec x2, int l2) const {
    double diff[NDIM];
    getBeadsSeparation(x1, l1, x2, l2, diff);

    double dist_sqrd = 0.0;

    for (int axis = 0; axis < NDIM; ++axis) {
        dist_sqrd += diff[axis] * diff[axis];
    }

    return dist_sqrd;
}

/* ---------------------------------------------------------------------- */

void BosonicExchangeBase::springForce(dVec& f) {
    if (bead_num == nbeads - 1) {
        springForceLastBead(f);
    } else if (bead_num == 0) {
        springForceFirstBead(f);
    } else {
        springForceInteriorBead(f);
    }
}

/**
 * Calculates the spring energy contribution of the current and the next time-slice,
 * provided that the connection is interior (i.e., not the spring energy between 1st and last bead).
 * This method cannot be called from the last time-slice.
 * 
 * @return Spring energy contribution of the current and the next time-slice
 */
double BosonicExchangeBase::interiorSpringEnergy() const {
    if (bead_num == nbeads - 1)
        throw std::runtime_error("interiorSpringEnergy() called for the last time-slice");

    double interior_spring_energy = 0.0;

    for (int ptcl_idx = 0; ptcl_idx < nbosons; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double diff = x(ptcl_idx, axis) - x_next(ptcl_idx, axis);
            interior_spring_energy += diff * diff;
        }
    }

    /// @todo Check if division by nbeads is necessary if i-Pi convention is used
    interior_spring_energy *= 0.5 * spring_constant;

    return interior_spring_energy;
}

/* ---------------------------------------------------------------------- */

void BosonicExchangeBase::springForceInteriorBead(dVec& f) const {
    for (int l = 0; l < nbosons; l++) {
        std::vector<double> sums(NDIM, 0.0);

        double diff_prev[NDIM];
        getBeadsSeparation(x, l, x_prev, l, diff_prev);

        for (int axis = 0; axis < NDIM; ++axis) {
            sums[axis] += diff_prev[axis];
        }

        double diff_next[3];
        getBeadsSeparation(x, l, x_next, l, diff_next);

        for (int axis = 0; axis < NDIM; ++axis) {
            sums[axis] += diff_next[axis];
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) = sums[axis] * spring_constant;
        }
    }
}
