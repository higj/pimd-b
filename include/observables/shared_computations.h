#pragma once

#include "contexts/velocity_context.h"
#include "contexts/bead_context.h"
#include "contexts/box_context.h"

class ObservableCache;
class ForceManager;

namespace SharedComputations {

    /**
     * @brief Returns the classical kinetic energy (internal units).
     * On the first call in a step it computes the value and caches it.
     * On subsequent calls it returns the cached result immediately.
     */
    double classicalKineticEnergy(
        const VelocityContext& vel_ctx,
        const BeadContext& bead_ctx,
        ObservableCache* cache   // nullptr => always compute (e.g. in tests)
    );

    /**
     * @brief Returns the raw external potential energy (pre /nbeads, internal units).
     * Result is cached under CacheKey::EXT_POT_RAW.
     */
    double extPotentialRaw(
        const VecArray& coord,
        const ForceManager& force_mgr,
        ObservableCache* cache
    );

    /**
     * @brief Returns the raw interaction potential energy (pre /nbeads, internal units).
     * Runs the O(N^2) pair loop on miss; result cached under CacheKey::INT_POT_RAW.
     */
    double intPotentialRaw(
        const VecArray& coord,
        const ForceManager& force_mgr,
        const BoxContext& box_ctx,
        const BeadContext& bead_ctx,
        ObservableCache* cache
    );

    /**
     * @brief Returns the raw virial sum -sum_i(r_i * F_i) (pre *0.5/nbeads).
     * Result is cached under CacheKey::VIRIAL_RAW.
     */
    double virialRaw(
        const VecArray& coord,
        const VecArray& forces,
        const BeadContext& bead_ctx,
        ObservableCache* cache
    );

    struct SuzukiChinComponents {
        double total_potential;
        double force_squared;
        double virial;
    };

    /**
     * @brief Computes (or retrieves from cache) the three raw intermediate values
     * shared by all Suzuki-Chin observables.
     *
     * Runs the O(N^2) pair loop - evaluating both V and gradV per pair - exactly once
     * per step. All three results are populated atomically; any observable that calls
     * this function after the first will get a full cache hit.
     */
    SuzukiChinComponents suzukiChinComponents(
        const VecArray& coord,
        const ForceManager& force_mgr,
        const BeadContext& bead_ctx,
        const BoxContext& box_ctx,
        ObservableCache* cache
    );

} // namespace SharedComputations