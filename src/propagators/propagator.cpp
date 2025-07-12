#include "propagators/propagator.h"
#include "core/system_state.h"
//#include "core/force_manager.h"
//#include "core/exchange_state.h"

Propagator::Propagator(
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<ForceManager>& force_mgr,
    const std::shared_ptr<ExchangeState>& exchange_state,
    const SpringContext& spring_ctx,
    double mass,
    double dt
) : m_state(state),
    m_force_mgr(force_mgr),
    m_exchange_state(exchange_state),
    m_spring_ctx(spring_ctx),
    m_natoms(state->getNumAtoms()),
    m_nbeads(state->getNumBeads()),
    m_this_bead(state->currentBead()),
    m_mass(mass),
    m_dt(dt)
{
}