#pragma once

#include "common.h"
#include "contexts/bosonic_exchange_context.h"

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
    /// TODO: Make context accept x_first_bead, x_last_bead, x_neighbor_bead as parameters, so that we can avoid assignFirstLast
    explicit BosonicExchangeBase(const BosonicExchangeContext& context);
    virtual ~BosonicExchangeBase() = default;
    ///BosonicExchangeBase(const BosonicExchangeBase&) = delete;
    ///BosonicExchangeBase& operator=(const BosonicExchangeBase&) = delete;

    void exteriorSpringForce(dVec& f);

    virtual void prepare() = 0;
    virtual double effectivePotential() = 0;
    virtual double primitiveEnergyEstimator() = 0;

    virtual double getDistinctProbability() = 0;
    virtual double getLongestProbability() = 0;

    virtual void printBosonicDebug() = 0;

protected:
    //void assignFirstLast(dVec& x_first_bead, dVec& x_last_bead) const;
    void getBeadsSeparation(const dVec& x1, int l1, const dVec& x2, int l2, double diff[NDIM]) const;
    [[nodiscard]] double getBeadsSeparationSquared(const dVec& x1, int l1, const dVec& x2, int l2) const;
    /// TODO: The interior calculation is slightly problematic, because it can involve either x_first_bead and x_last_bead
    /// TODO: Maybe move all interior spring calculations outside of this class?
    [[nodiscard]] double getInteriorSeparationSquared(int exterior_idx, int interior_idx) const;
    [[nodiscard]] double getExteriorSeparationSquared(int first_idx, int last_idx) const;

    // Pure virtual functions (must be implemented by derived classes)
    virtual void springForceFirstBead(dVec& f) = 0;
    virtual void springForceLastBead(dVec& f) = 0;

    BosonicExchangeContext m_context;
};
