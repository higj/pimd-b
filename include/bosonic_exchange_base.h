#pragma once

#include "common.h"

class BosonicExchangeBase {
public:
    BosonicExchangeBase(int nbosons, int np, int bead_num, double beta, double spring_constant,
        const dVec x, const dVec x_prev, const dVec x_next, bool pbc, double L);
    ~BosonicExchangeBase();

    virtual void updateCoordinates(const dVec new_x, const dVec new_x_prev, const dVec new_x_next);
    void spring_force(dVec& f);

    virtual double classical_potential();
    virtual double prim_estimator();

protected:
    void diff_two_beads(const dVec x1, int l1, const dVec x2, int l2, double diff[NDIM]);
    double distance_squared_two_beads(const dVec x1, int l1, const dVec x2, int l2);

    void spring_force_interior_bead(dVec& f);
    virtual void spring_force_first_bead(dVec& f);
    virtual void spring_force_last_bead(dVec& f);

    const int nbosons;
    const int np;
    const int bead_num;

    double spring_constant;
    double beta;
    dVec x;
    dVec x_prev;
    dVec x_next;

    bool pbc;
    double L; // Linear size of the system
};