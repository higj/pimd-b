#pragma once

#include "observables/observable.h"
#include "contexts/velocity_context.h"
#include "contexts/bead_context.h"

/**
 * @brief Observable for the classical kinetic energy of the ring polymers.
 *
 * Computes  K = 0.5 * sum_i( p_i^2 ) / mass  and records it as "cl_kinetic".
 *
 * The raw kinetic energy is stored in the ObservableCache under
 * CacheKey::CL_KINETIC_RAW so that TemperatureObservable, if active in the
 * same logging step, can reuse the result without repeating the momentum loop.
 */
class ClassicalKineticEnergyObservable final : public Observable {
public:
    ClassicalKineticEnergyObservable(
        const VelocityContext& vel_ctx,
        const BeadContext& bead_ctx,
        int   out_freq,
        const std::string& out_unit
    );

    void calculate() override;

private:
    VelocityContext m_vel_ctx;
    BeadContext     m_bead_ctx;
};