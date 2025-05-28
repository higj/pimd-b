#include <cassert>
#include <array>

#include "bosonic_exchange/bosonic_exchange_base.h"

BosonicExchangeBase::BosonicExchangeBase(const BosonicExchangeContext& context) : m_context(context) {
}

void BosonicExchangeBase::getExteriorBeadsSeparation(int first_idx, int last_idx, std::array<double, NDIM>& diff) const {
    first_idx = first_idx % m_context.nbosons;
    last_idx = last_idx % m_context.nbosons;

    for (int axis = 0; axis < NDIM; ++axis) {
        double dx = (*m_context.x_first_bead)(first_idx, axis) - (*m_context.x_last_bead)(last_idx, axis);
#if MINIM
        if (m_context.pbc)
            applyMinimumImage(dx, m_context.box_size);
#endif
        diff[axis] = dx;
    }
}

/*
void BosonicExchangeBase::getBeadsSeparation(const dVec& x1, int l1, const dVec& x2, int l2, double diff[NDIM]) const {
    l1 = l1 % m_context.nbosons;
    l2 = l2 % m_context.nbosons;

    for (int axis = 0; axis < NDIM; ++axis) {
        double dx = x2(l2, axis) - x1(l1, axis);
#if MINIM
        if (m_context.pbc)
            applyMinimumImage(dx, m_context.box_size);
#endif
        diff[axis] = dx;
    }
}
*/
/*
double BosonicExchangeBase::getBeadsSeparationSquared(const dVec& x1, int l1, const dVec& x2, int l2) const {
    double diff[NDIM];
    getBeadsSeparation(x1, l1, x2, l2, diff);

    double dist_sqrd = 0.0;

    for (int axis = 0; axis < NDIM; ++axis) {
        dist_sqrd += diff[axis] * diff[axis];
    }

    return dist_sqrd;
}
*/

double BosonicExchangeBase::getExteriorSeparationSquared(int first_idx, int last_idx) const {
    std::array<double, NDIM> diff;
    getExteriorBeadsSeparation(first_idx, last_idx, diff);

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
    if (m_context.this_bead == 0) {
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
/*void BosonicExchangeBase::assignFirstLast(dVec& x_first_bead, dVec& x_last_bead) const {
    if (m_context.this_bead == 0) {
        x_first_bead = *m_context.x;
        x_last_bead = *m_context.x_prev;
    } else {
        x_first_bead = *m_context.x_next;
        x_last_bead = *m_context.x;
    }
}
*/