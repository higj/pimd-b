#pragma once

#include "common.h"

class Simulation; // Forward declaration

/**
 * @class BosonicExchangeBase
 * @brief Abstract base class for bosonic exchange algorithms.
 *
 * Serves as a template for bosonic path integral molecular dynamics
 * algorithms. Basic functionality that is common to all bosonic algorithms,
 * such as updating the coordinates or evaluating the bead separation
 * is implemented here.
 */
class BosonicExchangeBase {
public:
    BosonicExchangeBase(const Simulation& _sim);
    virtual ~BosonicExchangeBase() = default;
    BosonicExchangeBase(const BosonicExchangeBase&) = delete;
    BosonicExchangeBase& operator=(const BosonicExchangeBase&) = delete;

    virtual void prepare() = 0;
    virtual void exteriorSpringForce(dVec& f) = 0;
    virtual double effectivePotential() = 0;
    virtual double primEstimator() = 0;

    virtual double getDistinctProbability() = 0;
    virtual double getLongestProbability() = 0;

    virtual void printBosonicDebug() = 0;

protected:
    virtual void assignFirstLast(dVec& x_first_bead, dVec& x_last_bead) const;
    void getBeadsSeparation(const dVec& x1, int l1, const dVec& x2, int l2, double diff[NDIM]) const;
    [[nodiscard]] double getBeadsSeparationSquared(const dVec& x1, int l1, const dVec& x2, int l2) const;

    const Simulation& sim; // Reference to the simulation object

    const int nbosons;
    const int nbeads;

    double spring_constant;
    double beta;

    const dVec& x;
    const dVec& x_prev;
    const dVec& x_next;

    virtual void assignIndirectionCoords();
    dVec indirection_x;
    dVec indirection_x_prev;
    dVec indirection_x_next;
    std::vector<int> indexes; 
};
