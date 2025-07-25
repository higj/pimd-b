#include "strategies/propagators/distinguishable_normal_modes_momenta_strategy.h"
#include "core/system_state.h"

void DistinguishableNormalModesMomentaStrategy::momentaExternalForces(
    const std::shared_ptr<SystemState>& state,
    const SpringContext& /* spring_ctx */,
    double dt
) {
#if !IPI_CONVENTION
    const int nbeads = state->getNumBeads();
#endif

    for (int ptcl_idx = 0; ptcl_idx < state->getNumAtoms(); ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
#if IPI_CONVENTION
            state->momenta(ptcl_idx, axis) += 0.5 * dt * state->physical_forces(ptcl_idx, axis);
#else
            state->momenta(ptcl_idx, axis) += 0.5 * dt * state->physical_forces(ptcl_idx, axis) / nbeads;
#endif
        }
    }
}
