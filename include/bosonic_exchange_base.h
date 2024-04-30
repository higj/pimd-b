#pragma once

#include "common.h"


/**
 * @class BosonicExchangeBase
 * @brief Abstract base class for bosonic exchange algorithms.
 * 
 * Serves as a template for bosonic path integral molecular dynamics
 * algorithms. Basic functionality that is common to all bosonic algorithms,
 * such as updating the coordinates or evaluating the interior spring
 * forces (which are unaffected by permutations) is implemented here.
 */
class BosonicExchangeBase {
public:
    BosonicExchangeBase(int nbosons_, int np_, int bead_num_, double beta_, double spring_constant_,
                        const dVec& x_, const dVec& x_prev_, const dVec& x_next_, bool pbc_, double size_);
    virtual ~BosonicExchangeBase() = default;
    BosonicExchangeBase(const BosonicExchangeBase&) = delete;
    BosonicExchangeBase& operator=(const BosonicExchangeBase&) = delete;

    void springForce(dVec& f);

    virtual void prepare() = 0;
    virtual double effectivePotential() = 0;
    virtual double primEstimator() = 0;

protected:
    void getBeadsSeparation(const dVec& x1, int l1, const dVec& x2, int l2, double diff[NDIM]) const;
    double getBeadsSeparationSquared(const dVec& x1, int l1, const dVec& x2, int l2) const;

    double interiorSpringEnergy() const;
    void springForceInteriorBead(dVec& f) const;

    // Pure virtual functions (must be implemented by derived classes)
    virtual void springForceFirstBead(dVec& f) = 0;
    virtual void springForceLastBead(dVec& f) = 0;

    const int nbosons;
    const int nbeads;
    const int bead_num;

    double spring_constant;
    double beta;

    const dVec& x;
    const dVec& x_prev;
    const dVec& x_next;

    bool pbc;
    double size; // Linear size of the system
};
