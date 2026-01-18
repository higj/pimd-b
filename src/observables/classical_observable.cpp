#include "observables/classical_observable.h"
#include "thermostats/thermostat.h"
#include "core/statistics_manager.h"
#include "strategies/observables/classical_spring_energy_strategy.h"

ClassicalObservable::ClassicalObservable(
    const std::shared_ptr<const dVec>& coord,
    const std::shared_ptr<const dVec>& prev_coord,
    const VelocityContext& vel_ctx,
    const ThermostatContext& thermostat_ctx,
    const BeadContext& bead_ctx,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    int out_freq,
    const std::string& out_unit
) : Observable("classical", out_freq, out_unit),
    m_coord(coord),
    m_prev_coord(prev_coord),
    m_spring_energy_strategy(StatisticsManager::getInstance().createClassicalSpringEnergyStrategy()),
    m_vel_ctx(vel_ctx),
    m_thermostat_ctx(thermostat_ctx),
    m_bead_ctx(bead_ctx),
    m_spring_ctx(spring_ctx),
    m_box_ctx(box_ctx)
{
    m_is_nose_hoover = thermostat_ctx.thermostat_type.find("nose_hoover") != std::string::npos;
    //const bool is_nose_hoover = thermostat_ctx.thermostat_type == "nose_hoover" || thermostat_ctx.thermostat_type == "nose_hoover_np" || thermostat_ctx.thermostat_type == "nose_hoover_np_dim";

    initialize({ "temperature", "cl_kinetic", "cl_spring" });

    if (m_is_nose_hoover)
        initializeLabel("nh_energy");
}

ClassicalObservable::~ClassicalObservable() = default;

void ClassicalObservable::calculate() {
    calculateKineticEnergy();
    calculateSpringEnergy();
    calculateThermostatEnergy();
}

void ClassicalObservable::calculateKineticEnergy() {
    double kinetic_energy = 0.0;

    const dVec& momenta = *m_vel_ctx.momenta;

    for (int ptcl_idx = 0; ptcl_idx < m_bead_ctx.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            kinetic_energy += momenta(ptcl_idx, axis) * momenta(ptcl_idx, axis);
        }
    }

    kinetic_energy *= 0.5 / m_vel_ctx.mass;

    quantities["cl_kinetic"] = Units::convertToUser("energy", m_out_unit, kinetic_energy);

    // Temperature is calculated according to Tolman's equipartition theorem as the average kinetic 
    // energy per degree of freedom. This might not be accurate for systems with constraints.
    // See [J. Chem. Theory Comput. 2019, 15, 1, 84-94.] for a discussion on the topic.

    /// @todo When zeroing the center of mass motion, the number of degrees of freedom must be reduced by NDIM.
    double dof = NDIM * m_bead_ctx.natoms * m_bead_ctx.nbeads;
    double temperature = 2.0 * kinetic_energy / (dof * Constants::kB);

    // In the i-Pi convention, the ring-polymer simulation is performed at a temperature that is P times higher
    // than the actual (quantum) temperature. Therefore, to ensure the quantum temperature is calculated correctly,
    // one must divide the classical temperature by the number of beads.
#if IPI_CONVENTION
    temperature /= m_bead_ctx.nbeads;
#endif

    /// @todo Allow conversion to different temperature units
    quantities["temperature"] = Units::convertToUser("temperature", "kelvin", temperature);
}

void ClassicalObservable::calculateSpringEnergy() {
    const double spring_energy = m_spring_energy_strategy->calculateSpringEnergy(
        *m_coord, 
        *m_prev_coord,
        m_spring_ctx, 
        m_box_ctx
    );

    quantities["cl_spring"] = Units::convertToUser("energy", m_out_unit, spring_energy);
}

void ClassicalObservable::calculateThermostatEnergy() {
    if (m_is_nose_hoover) {
        quantities["nh_energy"] = Units::convertToUser("energy", m_out_unit, m_thermostat_ctx.thermostat->getAdditionToH());
    }
}