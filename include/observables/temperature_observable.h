#pragma once

#include "observables/observable.h"
#include "contexts/velocity_context.h"
#include "contexts/bead_context.h"

/**
 * @brief Observable for the instantaneous classical temperature of the system.
 *
 * Temperature is derived from the classical kinetic energy via the equipartition
 * theorem:  T = 2K / (f * k_B),  where f = NDIM * natoms * nbeads.
 *
 * When TAU_CONVENTION is active the classical ring-polymer ensemble runs at
 * P times the physical temperature, so the result is divided by nbeads to
 * recover the quantum temperature.
 *
 * The underlying kinetic energy sum is retrieved from the ObservableCache
 * (CacheKey::CL_KINETIC_RAW) if ClKineticObservable has already run this step,
 * otherwise it is computed here and placed in the cache for any later consumer.
 */
class TemperatureObservable final : public Observable {
public:
    TemperatureObservable(
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