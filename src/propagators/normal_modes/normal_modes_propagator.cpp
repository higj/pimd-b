#include "propagators/normal_modes/normal_modes_propagator.h"
#include "core/force_field_manager.h"

#include <numbers>

NormalModesPropagator::NormalModesPropagator(const PropagatorContext& context, const std::shared_ptr<NormalModes>& normal_modes) : Propagator(context), m_normal_modes(normal_modes)
{
    // Frequencies
    m_freq = 2 * m_context.omega_p * sin(m_context.this_bead * std::numbers::pi / m_context.nbeads);
    m_c = cos(m_freq * m_context.dt);
    m_s = sin(m_freq * m_context.dt);
    m_omega = m_context.mass * m_freq;

    /*m_normal_modes = std::make_unique<NormalModes>(
        NormalModesContext{
            //.coord = std::make_shared<const dVec>(context.state->coord),
            //.momenta = std::make_shared<dVec>(context.state->momenta),
            .coord = std::shared_ptr<const dVec>(context.state, &context.state->coord),
            .momenta = std::shared_ptr<dVec>(context.state, &context.state->momenta),
            .natoms = context.natoms,
            .nbeads = context.nbeads,
            .this_bead = context.this_bead,
        }
    );
    */
}

void NormalModesPropagator::step() {
    // Propagate momenta under external forces
    momentaExternalForces();
    
    m_normal_modes->shareData();

    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            const int glob_idx = m_normal_modes->globIndexAtom(axis, ptcl_idx);
            // Cartesian-to-nm transformation
            const double coord_nm = m_normal_modes->coordCartesianToNormal(glob_idx);
            const double momentum_nm = m_normal_modes->momentumCartesianToNormal(glob_idx);
            // Time propagation
            if (m_freq == 0) {
                m_normal_modes->arr_coord_nm[glob_idx + m_context.this_bead] = coord_nm + m_context.dt / m_context.mass * momentum_nm;
                m_normal_modes->arr_momenta_nm[glob_idx + m_context.this_bead] = momentum_nm;
            } else {
                m_normal_modes->arr_coord_nm[glob_idx + m_context.this_bead] = m_c * coord_nm + m_s / m_omega * momentum_nm;
                m_normal_modes->arr_momenta_nm[glob_idx + m_context.this_bead] = (-1) * m_omega * m_s * coord_nm + m_c * momentum_nm;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = m_normal_modes->globIndexAtom(axis, ptcl_idx);
            // NM-to-Cartesian transformation
            const double coord_cartesian = m_normal_modes->coordNormalToCartesian(glob_idx);
            const double momentum_cartesian = m_normal_modes->momentumNormalToCartesian(glob_idx);
            m_normal_modes->arr_coord_cartesian[glob_idx + m_context.this_bead] = coord_cartesian;
            m_normal_modes->arr_momenta_cartesian[glob_idx + m_context.this_bead] = momentum_cartesian;
        }
    }

    // Update forces
    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = m_normal_modes->globIndexAtom(axis, ptcl_idx);
            m_context.state->coord(ptcl_idx, axis) = m_normal_modes->arr_coord_cartesian[glob_idx + m_context.this_bead];
            m_context.state->momenta(ptcl_idx, axis) = m_normal_modes->arr_momenta_cartesian[glob_idx + m_context.this_bead];
        }
    }

    /*sim.updateNeighboringCoordinates();
    sim.updateForces();
    sim.updatePhysicalForces(ext_forces);
    sim.updateSpringForces(spring_forces);*/

    // Remember to update the neighboring coordinates after every coordinate propagation
    m_context.state->updateNeighboringCoordinates();

    // Third step: forces are updated using the new positions
    m_context.force_mgr->updateSpringForces(*m_context.state, *m_context.exchange_state);
    m_context.force_mgr->updatePhysicalForces(*m_context.state);
    
    // Propagate momenta under external forces
    momentaExternalForces();
}

void NormalModesPropagator::momentaExternalForces() const
{
    if (!m_context.bosonic) {
        for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
#if IPI_CONVENTION
                m_context.state->momenta(ptcl_idx, axis) += 0.5 * m_context.dt * m_context.state->physical_forces(ptcl_idx, axis);
#else
                m_context.state->momenta(ptcl_idx, axis) += 0.5 * m_context.dt * m_context.state->physical_forces(ptcl_idx, axis) / m_context.nbeads;
#endif
    } else if (m_context.this_bead == 0 || m_context.this_bead == m_context.nbeads - 1) {
        for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                double inner_springs = -m_context.spring_constant * (2 * m_context.state->coord(ptcl_idx, axis) - m_context.state->prev_coord(ptcl_idx, axis) - m_context.state->
                    next_coord(ptcl_idx, axis));
#if IPI_CONVENTION
                m_context.state->momenta(ptcl_idx, axis) += 0.5 * m_context.dt * (m_context.state->physical_forces(ptcl_idx, axis) + m_context.state->spring_forces(ptcl_idx, axis) - inner_springs);
#else
                m_context.state->momenta(ptcl_idx, axis) += 0.5 * m_context.dt * (m_context.state->physical_forces(ptcl_idx, axis) / m_context.nbeads + m_context.state->spring_forces(ptcl_idx, axis) - inner_springs);
#endif
            }
        }
    } else {
        for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
#if IPI_CONVENTION
                m_context.state->momenta(ptcl_idx, axis) += 0.5 * m_context.dt * m_context.state->physical_forces(ptcl_idx, axis);
#else
                m_context.state->momenta(ptcl_idx, axis) += 0.5 * m_context.dt * m_context.state->physical_forces(ptcl_idx, axis) / m_context.nbeads;
#endif
    }
}
