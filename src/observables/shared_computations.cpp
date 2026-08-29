#include "observables/shared_computations.h"
#include "observables/observable_cache.h"
#include "observables/cache_keys.h"
#include "common.h"

double SharedComputations::classicalKineticEnergy(
    const VelocityContext& vel_ctx,
    const BeadContext& bead_ctx,
    ObservableCache* cache
) {
    if (cache) {
        if (const auto cached = cache->get(CacheKey::CL_KINETIC_RAW))
        {
            return *cached;
        }
    }

    double ke = 0.0;
    const VecArray& momenta = *vel_ctx.momenta;
    for (int ptcl_idx = 0; ptcl_idx < bead_ctx.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis)
            ke += momenta(ptcl_idx, axis) * momenta(ptcl_idx, axis);
    }
    ke *= 0.5 / vel_ctx.mass;

    if (cache) {
        cache->put(CacheKey::CL_KINETIC_RAW, ke);
    }

    return ke;
}