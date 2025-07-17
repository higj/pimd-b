#include "strategies/propagators/bosonic_normal_modes_momenta_strategy.h"
#include "core/system_state.h"

BosonicNormalModesMomentaStrategy::BosonicNormalModesMomentaStrategy(
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange
) : m_bosonic_exchange(bosonic_exchange) {
}

void BosonicNormalModesMomentaStrategy::momentaExternalForces(
    const std::shared_ptr<SystemState>& state,
    const SpringContext& spring_ctx,
    double dt
) {
    const int nbeads = state->getNumBeads();

    for (int ptcl_idx = 0; ptcl_idx < state->getNumAtoms(); ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double inner_springs = -spring_ctx.spring_constant * (2 * state->coord(ptcl_idx, axis) - state->prev_coord(ptcl_idx, axis) - state->next_coord(ptcl_idx, axis));
#if IPI_CONVENTION
            state->momenta(ptcl_idx, axis) += 0.5 * dt * (state->physical_forces(ptcl_idx, axis) + state->spring_forces(ptcl_idx, axis) - inner_springs);
#else
            state->momenta(ptcl_idx, axis) += 0.5 * dt * (state->physical_forces(ptcl_idx, axis) / nbeads + state->spring_forces(ptcl_idx, axis) - inner_springs);
#endif
        }
    }
}
