#include "observables/shared_computations.h"
#include "observables/observable_cache.h"
#include "observables/cache_keys.h"
#include "core/force_manager.h"  // for ext_potential, int_potential, cutoff
#include "common.h"

double SharedComputations::classicalKineticEnergy(
    const VelocityContext& vel_ctx,
    const BeadContext& bead_ctx,
    ObservableCache* cache
)
{
    if (cache)
    {
        if (const auto cached = cache->get(CacheKey::CL_KINETIC_RAW))
        {
            return *cached;
        }
    }

    double ke = 0.0;
    const VecArray& momenta = *vel_ctx.momenta;
    for (int ptcl_idx = 0; ptcl_idx < bead_ctx.natoms; ++ptcl_idx)
    {
        for (int axis = 0; axis < NDIM; ++axis)
            ke += momenta(ptcl_idx, axis) * momenta(ptcl_idx, axis);
    }
    ke *= 0.5 / vel_ctx.mass;

    if (cache)
    {
        cache->put(CacheKey::CL_KINETIC_RAW, ke);
    }

    return ke;
}

double SharedComputations::extPotentialRaw(
    const VecArray& coord,
    const ForceManager& force_mgr,
    ObservableCache* cache
)
{
    if (cache)
    {
        if (const auto cached = cache->get(CacheKey::EXT_POT_RAW))
            return *cached;
    }

    const double result = force_mgr.ext_potential->V(coord);

    if (cache)
    {
        cache->put(CacheKey::EXT_POT_RAW, result);
    }

    return result;
}

double SharedComputations::intPotentialRaw(
    const VecArray& coord,
    const ForceManager& force_mgr,
    const BoxContext& box_ctx,
    const BeadContext& bead_ctx,
    ObservableCache* cache
)
{
    if (cache)
    {
        if (const auto cached = cache->get(CacheKey::INT_POT_RAW))
            return *cached;
    }

    double result = 0.0;
    if (force_mgr.cutoff != 0.0)
    {
        for (int i = 0; i < bead_ctx.natoms; ++i)
        {
            for (int j = i + 1; j < bead_ctx.natoms; ++j)
            {
                Vec diff = coord.getSeparationArray(i, j);
                box_ctx.applyMinimumImageIfNeeded(diff);
                const double dist = norm(diff);
                if (dist < force_mgr.cutoff || force_mgr.cutoff < 0.0)
                    result += force_mgr.int_potential->V(diff);
            }
        }
    }

    if (cache)
    {
        cache->put(CacheKey::INT_POT_RAW, result);
    }

    return result;
}

double SharedComputations::virialRaw(
    const VecArray& coord,
    const VecArray& forces,
    const BeadContext& bead_ctx,
    ObservableCache* cache
)
{
    if (cache)
    {
        if (const auto cached = cache->get(CacheKey::VIRIAL_RAW))
        {
            return *cached;
        }
    }

    double result = 0.0;
    for (int ptcl_idx = 0; ptcl_idx < bead_ctx.natoms; ++ptcl_idx)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            result -= coord(ptcl_idx, axis) * forces(ptcl_idx, axis);
        }
    }

    if (cache)
    {
        cache->put(CacheKey::VIRIAL_RAW, result);
    }

    return result;
}

SharedComputations::SuzukiChinComponents SharedComputations::suzukiChinComponents(
    const VecArray& coord,
    const ForceManager& force_mgr,
    const BeadContext& bead_ctx,
    const BoxContext& box_ctx,
    ObservableCache* cache
)
{
    // Cache hit: all three values are always stored together, so checking one is enough.
    if (cache)
    {
        if (const auto cached_pot = cache->get(CacheKey::SC_TOTAL_POTENTIAL))
        {
            return {
                .total_potential = *cached_pot,
                .force_squared = *cache->get(CacheKey::SC_FORCE_SQUARED),
                .virial = *cache->get(CacheKey::SC_VIRIAL)
            };
        }
    }

    // External potential and its gradients
    double total_potential = force_mgr.ext_potential->V(coord);
    VecArray gradients(bead_ctx.natoms);
    gradients = force_mgr.ext_potential->gradV(coord);

    // External potential contribution to virial: sum_i coord_i * gradV_ext_i
    double virial = 0.0;
    if (!force_mgr.ext_potential->isFree())
    {
        for (int ptcl_idx = 0; ptcl_idx < bead_ctx.natoms; ++ptcl_idx)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                virial += coord(ptcl_idx, axis) * gradients(ptcl_idx, axis);
            }
        }
    }

    // Interaction potential pair loop: V, gradV, and virial contribution
    if (force_mgr.cutoff != 0.0 && !force_mgr.int_potential->isFree())
    {
        for (int i = 0; i < bead_ctx.natoms; ++i)
        {
            for (int j = i + 1; j < bead_ctx.natoms; ++j)
            {
                // TODO: apply minimum image (same as in GSFActionObservable)
                Vec diff = coord.getSeparationArray(i, j);
                if (const double dist = norm(diff); dist < force_mgr.cutoff || force_mgr.cutoff < 0.0)
                {
                    total_potential += force_mgr.int_potential->V(diff);
                    Vec int_grad = force_mgr.int_potential->gradV(diff);
                    box_ctx.applyMinimumImageIfNeeded(int_grad);

                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        gradients(i, axis) += int_grad[axis];
                        gradients(j, axis) -= int_grad[axis];
                        virial += coord(i, axis) * int_grad[axis]; // r_i * gradV(r_ij)
                        virial -= coord(j, axis) * int_grad[axis]; // r_j * gradV(r_ij)
                    }
                }
            }
        }
    }

    // Gradients squared
    double force_squared = 0.0;
    for (int ptcl_idx = 0; ptcl_idx < bead_ctx.natoms; ++ptcl_idx)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            force_squared += gradients(ptcl_idx, axis) * gradients(ptcl_idx, axis);
        }
    }

    if (cache)
    {
        cache->put(CacheKey::SC_TOTAL_POTENTIAL, total_potential);
        cache->put(CacheKey::SC_FORCE_SQUARED, force_squared);
        cache->put(CacheKey::SC_VIRIAL, virial);
    }

    return {
        .total_potential = total_potential,
        .force_squared = force_squared,
        .virial = virial
    };
}
