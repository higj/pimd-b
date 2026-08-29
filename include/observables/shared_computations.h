#pragma once

#include "contexts/velocity_context.h"
#include "contexts/bead_context.h"

class ObservableCache;

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

} // namespace SharedComputations