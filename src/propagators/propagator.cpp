#include "propagators/propagator.h"
#include "core/system_state.h"

Propagator::Propagator(
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<ForceManager>& force_mgr,
    const BoxContext& box_ctx,
    const SpringContext& spring_ctx,
    double mass,
    double dt
) : m_state(state),
    m_force_mgr(force_mgr),
    m_box_ctx(box_ctx),
    m_spring_ctx(spring_ctx),
    m_natoms(state->getNumAtoms()),
    m_nbeads(state->getNumBeads()),
    m_this_bead(state->currentBead()),
    m_mass(mass),
    m_dt(dt)
{
}