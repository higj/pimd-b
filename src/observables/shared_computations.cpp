#include "observables/shared_computations.h"
#include "observables/observable_cache.h"
#include "observables/cache_keys.h"
#include "core/force_manager.h"  // for ext_potential, int_potential, cutoff
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

double SharedComputations::extPotentialRaw(
    const VecArray& coord,
    const ForceManager& force_mgr,
    ObservableCache* cache
) {
    if (cache) {
        if (const auto cached = cache->get(CacheKey::EXT_POT_RAW))
            return *cached;
    }
    const double result = force_mgr.ext_potential->V(coord);
    if (cache) cache->put(CacheKey::EXT_POT_RAW, result);
    return result;
}

double SharedComputations::intPotentialRaw(
    const VecArray& coord,
    const ForceManager& force_mgr,
    const BoxContext& box_ctx,
    const BeadContext& bead_ctx,
    ObservableCache* cache
) {
    if (cache) {
        if (const auto cached = cache->get(CacheKey::INT_POT_RAW))
            return *cached;
    }
    double result = 0.0;
    if (force_mgr.cutoff != 0.0) {
        for (int i = 0; i < bead_ctx.natoms; ++i) {
            for (int j = i + 1; j < bead_ctx.natoms; ++j) {
                Vec diff = coord.getSeparationArray(i, j);
                box_ctx.applyMinimumImageIfNeeded(diff);
                const double dist = norm(diff);
                if (dist < force_mgr.cutoff || force_mgr.cutoff < 0.0)
                    result += force_mgr.int_potential->V(diff);
            }
        }
    }
    if (cache) cache->put(CacheKey::INT_POT_RAW, result);
    return result;
}

double SharedComputations::virialRaw(
    const VecArray& coord,
    const VecArray& forces,
    const BeadContext& bead_ctx,
    ObservableCache* cache
) {
    if (cache) {
        if (const auto cached = cache->get(CacheKey::VIRIAL_RAW))
            return *cached;
    }
    double result = 0.0;
    for (int i = 0; i < bead_ctx.natoms; ++i)
        for (int a = 0; a < NDIM; ++a)
            result -= coord(i, a) * forces(i, a);
    if (cache) cache->put(CacheKey::VIRIAL_RAW, result);
    return result;
}