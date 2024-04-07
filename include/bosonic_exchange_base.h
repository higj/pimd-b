#pragma once

#include "common.h"

class BosonicExchangeBase {
public:
    BosonicExchangeBase(int nbosons, int np, int bead_num, double beta, double spring_constant,
        const dVec x, const dVec x_prev, const dVec x_next, bool pbc, double L);
    ~BosonicExchangeBase();

    virtual void updateCoordinates(const dVec new_x, const dVec new_x_prev, const dVec new_x_next);
    void springForce(dVec& f);
    virtual double effectivePotential();

    virtual double primEstimator();

protected:
    void diff_two_beads(const dVec x1, int l1, const dVec x2, int l2, double diff[NDIM]);
    double distance_squared_two_beads(const dVec x1, int l1, const dVec x2, int l2);

    double interiorSpringEnergy();
    void springForceInteriorBead(dVec& f);
    virtual void springForceFirstBead(dVec& f);
    virtual void springForceLastBead(dVec& f);

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