#include "strategies/forces/distinguishable_spring_force_strategy.h"

void DistinguishableSpringForceStrategy::updateSpringForces(SystemState& state, const SpringContext& spring_ctx, const BoxContext& box_ctx) {
    for (int ptcl_idx = 0; ptcl_idx < state.getNumAtoms(); ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double diff_prev = state.prev_coord(ptcl_idx, axis) - state.coord(ptcl_idx, axis);
            double diff_next = state.next_coord(ptcl_idx, axis) - state.coord(ptcl_idx, axis);

            box_ctx.applyMinimumImageIfNeeded(diff_prev);
            box_ctx.applyMinimumImageIfNeeded(diff_next);

            state.spring_forces(ptcl_idx, axis) = spring_ctx.spring_constant * (diff_prev + diff_next);
        }
    }
}