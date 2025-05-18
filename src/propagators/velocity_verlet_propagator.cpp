#include "propagators/velocity_verlet_propagator.h"

#include "core/system_state.h"
#include "core/force_field_manager.h"

VelocityVerletPropagator::VelocityVerletPropagator(const PropagatorContext& context) : Propagator(context) {
}

void VelocityVerletPropagator::step() {
    // First step: momenta are propagated by half a step ("B" step)
    momentStep();

    // Second step: positions are propagated using the new momenta ("A" step)
    coordsStep();

    // Remember to update the neighboring coordinates after every coordinate propagation
    m_context.state->updateNeighboringCoordinates();

    // Third step: forces are updated using the new positions
    /// TODO: Passing state and exchange_state every time we need to update forces is not ideal (delegate this to the force manager)?
    /// TODO: It could be simpler to have m_context.force_mgr->updateForces();
    m_context.force_mgr->updateSpringForces(*m_context.state, *m_context.exchange_state);
    m_context.force_mgr->updatePhysicalForces(*m_context.state);

    // Fourth step: momenta are propagated once more ("B" step)
    momentStep();
}

void Propagator::momentStep() const
{
    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            m_context.state->momenta(ptcl_idx, axis) += 0.5 * m_context.dt * m_context.state->getTotalForce(ptcl_idx, axis);
        }
    }
}

void Propagator::coordsStep() const
{
    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            m_context.state->coord(ptcl_idx, axis) += m_context.dt * m_context.state->momenta(ptcl_idx, axis) / m_context.mass;
        }
    }
}
