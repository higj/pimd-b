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
    if (!m_state)
    {
        throw std::invalid_argument("SystemState pointer is null.");
    }

    if (!m_force_mgr)
    {
        throw std::invalid_argument("ForceManager pointer is null.");
    }

    if (m_natoms <= 0)
    {
        throw std::invalid_argument("Number of atoms must be positive.");
    }

    if (m_nbeads <= 0)
    {
        throw std::invalid_argument("Number of beads must be positive.");
    }

    if (m_mass <= 0.0)
    {
        throw std::invalid_argument("Mass must be positive.");
    }

    if (m_dt <= 0.0)
    {
        throw std::invalid_argument("Time step (dt) must be positive.");
    }
}