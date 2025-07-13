#pragma once

#include "common.h"
#include "contexts/box_context.h"
#include "contexts/thermal_context.h"
#include "contexts/spring_context.h"
#include "contexts/bead_context.h"

#include <memory>

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
    explicit BosonicExchangeBase(
        const std::shared_ptr<const dVec>& coord_first_bead,
        const std::shared_ptr<const dVec>& coord_last_bead,
        const ThermalContext& thermal_ctx,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        const BeadContext& bead_ctx
    );
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
    /**
     * Calculates the vector distance between two beads of an exterior spring (first minus last).
     *
     * @param first_idx Particle index of the bead at the first imaginary time-slice.
     * @param last_idx Particle index of the bead at the last imaginary time-slice.
     * @param diff Output array to store the distance vector.
     */
    void getExteriorBeadsSeparation(int first_idx, int last_idx, std::array<double, NDIM>& diff) const;

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

    std::shared_ptr<const dVec> m_coord_first_bead;
    std::shared_ptr<const dVec> m_coord_last_bead;

    ThermalContext m_thermal_ctx;
    SpringContext m_spring_ctx;
    BoxContext m_box_ctx;
    BeadContext m_bead_ctx;
};
