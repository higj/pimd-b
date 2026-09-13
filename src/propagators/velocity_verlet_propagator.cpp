#include "propagators/velocity_verlet_propagator.h"

#include "core/system_state.h"
#include "core/force_manager.h"

VelocityVerletPropagator::VelocityVerletPropagator(
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<ForceManager>& force_mgr,
    const BoxContext& box_ctx,
    const SpringContext& spring_ctx,
    double mass,
    double dt
) : Propagator(state, force_mgr, box_ctx, spring_ctx, mass, dt)
{
}

void VelocityVerletPropagator::step()
{
    // First step: momenta are propagated by half a step ("B" step)
    momentStep();

    // Second step: positions are propagated using the new momenta ("A" step)
    coordsStep();

    // Third step: forces are updated using the new positions
    m_force_mgr->updateForces(*m_state, m_spring_ctx, m_box_ctx);

    // Fourth step: momenta are propagated once more ("B" step)
    momentStep();
}

void VelocityVerletPropagator::momentStep() const
{
    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            m_state->momenta(ptcl_idx, axis) += 0.5 * m_dt * m_state->getTotalForce(ptcl_idx, axis);
        }
    }
}

void VelocityVerletPropagator::coordsStep() const
{
    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            m_state->coord(ptcl_idx, axis) += m_dt * m_state->momenta(ptcl_idx, axis) / m_mass;
        }
    }

    // Apply periodic boundary conditions after coordinate update
    m_box_ctx.applyMinimumImageIfNeeded(m_state->coord);

    // Remember to update the neighboring coordinates after every coordinate propagation
    m_state->updateNeighboringCoordinates();
}
