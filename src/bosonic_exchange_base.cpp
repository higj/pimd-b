#include <cassert>

#include "bosonic_exchange_base.h"
#include "simulation.h"

BosonicExchangeBase::BosonicExchangeBase(const Simulation& _sim) :
    sim(_sim),
    nbosons(_sim.natoms),
    nbeads(_sim.nbeads),
    spring_constant(_sim.spring_constant),
    beta(_sim.beta),
    x(_sim.coord),
    x_prev(_sim.prev_coord),
    x_next(_sim.next_coord) {
    assert(_sim.is_bosonic_bead);
#if IPI_CONVENTION
    beta /= _sim.nbeads;
#endif
}

/**
 * Calculates the vector distance between two beads (second minus first).
 *
 * @param x1 Coordinates of the first bead.
 * @param l1 Particle index of the first bead.
 * @param x2 Coordinates of the second bead.
 * @param l2 Particle index of the second bead.
 * @param[out] diff Vector distance between the two beads.
 */
void BosonicExchangeBase::getBeadsSeparation(const dVec& x1, int l1, const dVec& x2, int l2, double diff[NDIM]) const {
    l1 = l1 % nbosons;
    l2 = l2 % nbosons;

    for (int axis = 0; axis < NDIM; ++axis) {
        double dx = x2(l2, axis) - x1(l1, axis);
#if MINIM
        if (sim.pbc)
            applyMinimumImage(dx, sim.size);
#endif
        diff[axis] = dx;
    }
}

/**
 * Calculates the distance squared between two beads.
 * 
 * @param x1 Coordinates of the first bead.
 * @param l1 Particle index of the first bead.
 * @param x2 Coordinates of the second bead.
 * @param l2 Particle index of the second bead.
 * @return Distance squared between the two beads.
 */
double BosonicExchangeBase::getBeadsSeparationSquared(const dVec& x1, int l1, const dVec& x2, int l2) const {
    double diff[NDIM];
    getBeadsSeparation(x1, l1, x2, l2, diff);

    double dist_sqrd = 0.0;

    for (int axis = 0; axis < NDIM; ++axis) {
        dist_sqrd += diff[axis] * diff[axis];
    }

    return dist_sqrd;
}

/**
 * Calculates the bosonic force exerted on the beads
 * at the first and the last time-slices.
 *
 * @param[out] f Vector to store the forces.
 */
void BosonicExchangeBase::exteriorSpringForce(dVec& f) {
    if (sim.this_bead == 0) {
        springForceFirstBead(f);
    } else {
        springForceLastBead(f);
    }
}

/**
 * Determines the coordinates of the first and the last bead based on the current time-slice.
 *
 * @param[out] x_first_bead The coordinates of the particles in the first time-slice.
 * @param[out] x_last_bead The coordinates of the particles in the last time-slice.
 */
void BosonicExchangeBase::assignFirstLast(dVec& x_first_bead, dVec& x_last_bead) const {
    if (sim.this_bead == 0) {
        x_first_bead = x;
        x_last_bead = x_prev;
    } else {
        x_first_bead = x_next;
        x_last_bead = x;
    }
}
