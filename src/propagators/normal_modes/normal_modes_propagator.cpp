#include "propagators/normal_modes/normal_modes_propagator.h"
#include "core/force_manager.h"

#include <numbers>

NormalModesPropagator::NormalModesPropagator(
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<ForceManager>& force_mgr,
    const std::shared_ptr<ExchangeState>& exchange_state,
    const SpringContext& spring_ctx,
    double mass,
    double dt,
    const std::shared_ptr<NormalModes>& normal_modes
) : Propagator(state, force_mgr, exchange_state, spring_ctx, mass, dt), m_normal_modes(normal_modes)
{
    // Frequencies
    m_freq = 2 * spring_ctx.omega_p * sin(m_this_bead * std::numbers::pi / m_nbeads);
    m_c = cos(m_freq * m_dt);
    m_s = sin(m_freq * m_dt);
    m_omega = mass * m_freq;
}

void NormalModesPropagator::step() {
    // Propagate momenta under external forces
    momentaExternalForces();
    
    m_normal_modes->shareData();

    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            const int glob_idx = m_normal_modes->globIndexAtom(axis, ptcl_idx);
            // Cartesian-to-nm transformation
            const double coord_nm = m_normal_modes->coordCartesianToNormal(glob_idx);
            const double momentum_nm = m_normal_modes->momentumCartesianToNormal(glob_idx);
            // Time propagation
            if (m_freq == 0) {
                m_normal_modes->arr_coord_nm[glob_idx + m_this_bead] = coord_nm + m_dt / m_mass * momentum_nm;
                m_normal_modes->arr_momenta_nm[glob_idx + m_this_bead] = momentum_nm;
            } else {
                m_normal_modes->arr_coord_nm[glob_idx + m_this_bead] = m_c * coord_nm + m_s / m_omega * momentum_nm;
                m_normal_modes->arr_momenta_nm[glob_idx + m_this_bead] = (-1) * m_omega * m_s * coord_nm + m_c * momentum_nm;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            const int glob_idx = m_normal_modes->globIndexAtom(axis, ptcl_idx);
            // NM-to-Cartesian transformation
            const double coord_cartesian = m_normal_modes->coordNormalToCartesian(glob_idx);
            const double momentum_cartesian = m_normal_modes->momentumNormalToCartesian(glob_idx);
            m_normal_modes->arr_coord_cartesian[glob_idx + m_this_bead] = coord_cartesian;
            m_normal_modes->arr_momenta_cartesian[glob_idx + m_this_bead] = momentum_cartesian;
        }
    }

    // Update coordinates and momenta
    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = m_normal_modes->globIndexAtom(axis, ptcl_idx);
            m_state->coord(ptcl_idx, axis) = m_normal_modes->arr_coord_cartesian[glob_idx + m_this_bead];
            m_state->momenta(ptcl_idx, axis) = m_normal_modes->arr_momenta_cartesian[glob_idx + m_this_bead];
        }
    }

    // Remember to update the neighboring coordinates after every coordinate propagation
    m_state->updateNeighboringCoordinates();

    // Third step: forces are updated using the new positions
    m_force_mgr->updateForces(*m_state, *m_exchange_state);
    
    // Propagate momenta under external forces
    momentaExternalForces();
}

void NormalModesPropagator::momentaExternalForces() const
{
    if (!m_exchange_state->is_bosonic) {
        for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
#if IPI_CONVENTION
                m_state->momenta(ptcl_idx, axis) += 0.5 * m_dt * m_state->physical_forces(ptcl_idx, axis);
#else
                m_state->momenta(ptcl_idx, axis) += 0.5 * m_dt * m_state->physical_forces(ptcl_idx, axis) / m_nbeads;
#endif
    } else if (m_this_bead == 0 || m_this_bead == m_nbeads - 1) {
        for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                double inner_springs = -m_spring_ctx.spring_constant * (2 * m_state->coord(ptcl_idx, axis) - m_state->prev_coord(ptcl_idx, axis) - m_state->
                    next_coord(ptcl_idx, axis));
#if IPI_CONVENTION
                m_state->momenta(ptcl_idx, axis) += 0.5 * m_dt * (m_state->physical_forces(ptcl_idx, axis) + m_state->spring_forces(ptcl_idx, axis) - inner_springs);
#else
                m_state->momenta(ptcl_idx, axis) += 0.5 * m_dt * (m_state->physical_forces(ptcl_idx, axis) / m_nbeads + m_state->spring_forces(ptcl_idx, axis) - inner_springs);
#endif
            }
        }
    } else {
        for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
#if IPI_CONVENTION
                m_state->momenta(ptcl_idx, axis) += 0.5 * m_dt * m_state->physical_forces(ptcl_idx, axis);
#else
                m_state->momenta(ptcl_idx, axis) += 0.5 * m_dt * m_state->physical_forces(ptcl_idx, axis) / m_nbeads;
#endif
    }
}
