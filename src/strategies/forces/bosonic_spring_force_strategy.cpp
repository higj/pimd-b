#include "strategies/forces/bosonic_spring_force_strategy.h"
#include "bosonic_exchange/bosonic_exchange_base.h"

BosonicSpringForceStrategy::BosonicSpringForceStrategy(
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange, 
    const BeadContext& bead_ctx
) : m_bosonic_exchange(bosonic_exchange),
    m_bead_ctx(bead_ctx),
    m_is_exterior_bead(bead_ctx.this_bead == 0 || bead_ctx.this_bead == bead_ctx.nbeads - 1)
{
}

void BosonicSpringForceStrategy::updateSpringForces(SystemState& state, const SpringContext& spring_ctx, const BoxContext& box_ctx) {
    if (m_is_exterior_bead) {
        // If the simulation is bosonic and the current bead is either 1 or P, we calculate
        // the exterior spring forces in the appropriate bosonic class.
        updateBosonicForces(state, spring_ctx, box_ctx);
    } else {
        // If particles are distinguishable, or if the current bead is an interior bead,
        // the force is calculated based on the standard expression for distinguishable particles.
        DistinguishableSpringForceStrategy::updateSpringForces(state, spring_ctx, box_ctx);
    }
}

void BosonicSpringForceStrategy::updateBosonicForces(SystemState& state, const SpringContext& spring_ctx, const BoxContext& box_ctx)
{
    m_bosonic_exchange->prepare();
    m_bosonic_exchange->exteriorSpringForce(state.spring_forces);

    if (m_bead_ctx.nbeads == 1) {
        // For bosonic simulations with a single imaginary time slice (P=1),
        // no inter-slice springs exist - only permutation-related springs within the same slice.
        // These should already have been handled by exteriorSpringForce, so we can exit early.
        // In the distinguishable case, P=1 implies no springs at all, as only diagonal
        // elements of the density matrix are relevant.
        return;
    }

    // Bosonic exchange class only calculates the contribution due to the exterior spring.
    // However, beads 1 and P are also affected by the interior spring (due to beads 2 and P-1, respectively).
    if (m_bead_ctx.this_bead == 0) {
        for (int ptcl_idx = 0; ptcl_idx < m_bead_ctx.natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                double diff_next = state.next_coord(ptcl_idx, axis) - state.coord(ptcl_idx, axis);
                box_ctx.applyMinimumImageIfNeeded(diff_next);
                state.spring_forces(ptcl_idx, axis) += spring_ctx.spring_constant * diff_next;
            }
        }

        return;
    }

    // The following is the spring force contribution acting on the last bead, due to the spring between P-1 and P
    for (int ptcl_idx = 0; ptcl_idx < m_bead_ctx.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double diff_prev = state.prev_coord(ptcl_idx, axis) - state.coord(ptcl_idx, axis);
            box_ctx.applyMinimumImageIfNeeded(diff_prev);
            state.spring_forces(ptcl_idx, axis) += spring_ctx.spring_constant * diff_prev;
        }
    }
}
