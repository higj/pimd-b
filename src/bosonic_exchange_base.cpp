#include "bosonic_exchange_base.h"

BosonicExchangeBase::BosonicExchangeBase(int nbosons, int np, int bead_num, double beta, double spring_constant,
    const dVec x, const dVec x_prev, const dVec x_next, bool pbc, double L) : 
    nbosons(nbosons),
    np(np),
    bead_num(bead_num),
    beta(beta),
    spring_constant(spring_constant),
    x(x),
    x_prev(x_prev),
    x_next(x_next),
    pbc(pbc),
    L(L) {
#if IPI_CONVENTION
    this->beta = beta / np;
#endif
}

/* ---------------------------------------------------------------------- */

BosonicExchangeBase::~BosonicExchangeBase() {
}

/* ---------------------------------------------------------------------- */

void BosonicExchangeBase::updateCoordinates(const dVec new_x, const dVec new_x_prev, const dVec new_x_next) {
    x = new_x;
    x_prev = new_x_prev;
    x_next = new_x_next;
}

/* ---------------------------------------------------------------------- */

void BosonicExchangeBase::diff_two_beads(const dVec x1, int l1, const dVec x2, int l2, double diff[NDIM]) {
    l1 = l1 % nbosons;
    l2 = l2 % nbosons;

    for (int axis = 0; axis < NDIM; ++axis) {
        double dx = x2(l2, axis) - x1(l1, axis);
#if MINIM
        if (pbc)
            applyMinimumImage(dx, L);
#endif
        diff[axis] = dx;
    }
}

/* ---------------------------------------------------------------------- */

double BosonicExchangeBase::distance_squared_two_beads(const dVec x1, int l1, const dVec x2, int l2)
{
    double diff[NDIM];
    diff_two_beads(x1, l1, x2, l2, diff);

    double dist_sqrd = 0.0;

    for (int axis = 0; axis < NDIM; ++axis) {
        dist_sqrd += diff[axis] * diff[axis];
    }

    return dist_sqrd;
}

/* ---------------------------------------------------------------------- */

void BosonicExchangeBase::spring_force(dVec& f) {
    if (bead_num == np - 1) {
        spring_force_last_bead(f);
    }
    else if (bead_num == 0) {
        spring_force_first_bead(f);
    }
    else {
        spring_force_interior_bead(f);
    }
}

/* ---------------------------------------------------------------------- */

void BosonicExchangeBase::spring_force_interior_bead(dVec& f)
{
    for (int l = 0; l < nbosons; l++) {
        std::vector<double> sums(NDIM, 0.0);

        double diff_prev[NDIM];
        diff_two_beads(x, l, x_prev, l, diff_prev);

        for (int axis = 0; axis < NDIM; ++axis) {
            sums[axis] += diff_prev[axis];
        }

        double diff_next[3];
        diff_two_beads(x, l, x_next, l, diff_next);

        for (int axis = 0; axis < NDIM; ++axis) {
            sums[axis] += diff_next[axis];
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) = sums[axis] * spring_constant;
        }
    }
}

/* ---------------------------------------------------------------------- */

void BosonicExchangeBase::spring_force_first_bead(dVec& f)
{
}

/* ---------------------------------------------------------------------- */

void BosonicExchangeBase::spring_force_last_bead(dVec& f)
{
}

/* ---------------------------------------------------------------------- */

double BosonicExchangeBase::classical_potential()
{
    return 0.0;
}

/* ---------------------------------------------------------------------- */

double BosonicExchangeBase::prim_estimator()
{
	return 0.0;
}
