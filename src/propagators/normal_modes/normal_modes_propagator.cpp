#include "propagators/normal_modes/normal_modes_propagator.h"
#include "core/force_manager.h"
#include "core/statistics_manager.h"

#include <numbers>

NormalModesPropagator::NormalModesPropagator(
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<ForceManager>& force_mgr,
    const BoxContext& box_ctx,
    const SpringContext& spring_ctx,
    double mass,
    double dt,
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
    const BeadContext& bead_ctx,
    bool bosonic,
    const std::shared_ptr<NormalModes>& normal_modes
) : Propagator(state, force_mgr, box_ctx, spring_ctx, mass, dt),
    m_normal_modes(normal_modes),
    m_nm_momenta_strategy(StatisticsManager::createNormalModesMomentaStrategy(bosonic_exchange, bead_ctx, bosonic))
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
    m_force_mgr->updateForces(*m_state, m_spring_ctx, m_box_ctx);
    
    // Propagate momenta under external forces
    momentaExternalForces();
}

void NormalModesPropagator::momentaExternalForces() const
{
    m_nm_momenta_strategy->momentaExternalForces(m_state, m_spring_ctx, m_dt);
}
