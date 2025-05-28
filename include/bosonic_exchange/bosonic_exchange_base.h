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
    /**
     * Calculates the vector distance between two beads of an exterior spring (first minus last).
     *
     * @param first_idx Particle index of the bead at the first imaginary time-slice.
     * @param last_idx Particle index of the bead at the last imaginary time-slice.
     * @param diff Output array to store the distance vector.
     */
    void getExteriorBeadsSeparation(int first_idx, int last_idx, std::array<double, NDIM>& diff) const;
    void getBeadsSeparation(const dVec& x1, int l1, const dVec& x2, int l2, double diff[NDIM]) const;

    /**
     * Calculates the distance squared between two beads.
     *
     * @param x1 Coordinates of the first bead.
     * @param l1 Particle index of the first bead.
     * @param x2 Coordinates of the second bead.
     * @param l2 Particle index of the second bead.
     * @return Distance squared between the two beads.
     */
    //[[nodiscard]] double getBeadsSeparationSquared(const dVec& x1, int l1, const dVec& x2, int l2) const;

    /**
     * Calculates the distance squared between two beads of an exterior spring.
     *
     * @param first_idx Particle index of the bead at the first imaginary time-slice.
     * @param last_idx Particle index of the bead at the last imaginary time-slice.
     * @return Distance squared between the two beads.
     */
    [[nodiscard]] double getExteriorSeparationSquared(int first_idx, int last_idx) const;

    // Pure virtual functions (must be implemented by derived classes)
    virtual void springForceFirstBead(dVec& f) = 0;
    virtual void springForceLastBead(dVec& f) = 0;

    BosonicExchangeContext m_context;
};
