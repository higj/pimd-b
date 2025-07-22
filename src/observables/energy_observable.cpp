#include "observables/energy_observable.h"
#include "core/force_manager.h"
#include "core/statistics_manager.h"
#include "strategies/observables/primitive_kinetic_energy_strategy.h"

EnergyObservable::EnergyObservable(
        const std::shared_ptr<const dVec>& coord,
        const std::shared_ptr<const dVec>& prev_coord,
        const std::shared_ptr<const ForceManager>& force_mgr,
        const BeadContext& bead_ctx,
        const ThermalContext& thermal_ctx,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        int out_freq,
        const std::string& out_unit
    ) : Observable(out_freq, out_unit),
    m_prim_ke_strategy(StatisticsManager::getInstance().createPrimitiveKineticEnergyStrategy()),
    m_coord_this(coord),
    m_coord_prev(prev_coord),
    m_force_mgr(force_mgr),
    m_bead_ctx(bead_ctx),
    m_thermal_ctx(thermal_ctx),
    m_spring_ctx(spring_ctx),
    m_box_ctx(box_ctx),
    m_is_ext_free(force_mgr->ext_potential->isFree()),
    m_is_int_free(force_mgr->int_potential->isFree())
{
    initialize({"kinetic"});

    if (!m_is_ext_free || !m_is_int_free)
    {
        initialize({"potential", "virial"});

        if (!m_is_ext_free && !m_is_int_free)
        {
            initialize({"ext_pot", "int_pot"});
        }
    }
}

EnergyObservable::~EnergyObservable() = default;

void EnergyObservable::calculate()
{
    calculateKinetic();
    calculatePotential();
}

void EnergyObservable::calculateKinetic()
{
    // First, add the constant factor of d*N*P/(2*beta) to the kinetic energy (per bead)
    quantities["kinetic"] = 0.5 * NDIM * m_bead_ctx.natoms / m_thermal_ctx.beta;
    quantities["kinetic"] += m_prim_ke_strategy->calculateSpringContribution(
        *m_coord_this, 
        *m_coord_prev, 
        m_spring_ctx, 
        m_box_ctx, 
        m_bead_ctx
    );

    quantities["kinetic"] = Units::convertToUser("energy", m_out_unit, quantities["kinetic"]);
}

void EnergyObservable::calculatePotential()
{
    double potential = 0.0; // Total potential energy
    double virial = 0.0; // Virial kinetic energy
    double int_pot = 0.0; // Potential energy due to interactions
    double ext_pot = 0.0; // Potential energy due to external field

    const auto& coord = *m_coord_this;

    if (!m_is_ext_free)
    {
        ext_pot = m_force_mgr->ext_potential->V(coord);
        potential += ext_pot;

        /// TODO: We already have physical_forces in ForceManager, so we can use it instead of calculating it again?
        dVec physical_forces(m_bead_ctx.natoms);
        physical_forces = (-1.0) * m_force_mgr->ext_potential->gradV(coord);

        for (int ptcl_idx = 0; ptcl_idx < m_bead_ctx.natoms; ++ptcl_idx)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                virial -= coord(ptcl_idx, axis) * physical_forces(ptcl_idx, axis);
            }
        }
    }

    if (m_force_mgr->cutoff != 0.0)
    {
        for (int ptcl_one = 0; ptcl_one < m_bead_ctx.natoms; ++ptcl_one)
        {
            for (int ptcl_two = ptcl_one + 1; ptcl_two < m_bead_ctx.natoms; ++ptcl_two)
            {
                dVec diff = coord.getSeparation(ptcl_one, ptcl_two); // Vectorial distance
                m_box_ctx.applyMinimumImageIfNeeded(diff);

                if (const double distance = diff.norm(); distance < m_force_mgr->cutoff || m_force_mgr->cutoff < 0.0)
                {
                    dVec force_on_one = (-1.0) * m_force_mgr->int_potential->gradV(diff);

                    double int_pot_val = m_force_mgr->int_potential->V(diff);
                    potential += int_pot_val;
                    int_pot += int_pot_val;

                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        virial -= coord(ptcl_one, axis) * force_on_one(0, axis);
                    }
                }
            }
        }
    }

    if (!m_is_ext_free && !m_is_int_free)
    {
        ext_pot /= m_bead_ctx.nbeads;
        int_pot /= m_bead_ctx.nbeads;

        quantities["ext_pot"] = Units::convertToUser("energy", m_out_unit, ext_pot);
        quantities["int_pot"] = Units::convertToUser("energy", m_out_unit, int_pot);
    }

    if (!m_is_ext_free || !m_is_int_free)
    {
        potential /= m_bead_ctx.nbeads;
        virial *= 0.5 / m_bead_ctx.nbeads;

        quantities["potential"] = Units::convertToUser("energy", m_out_unit, potential);
        quantities["virial"] = Units::convertToUser("energy", m_out_unit, virial);
    }
}
