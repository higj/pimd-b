#include "observables/classical_observable.h"
#include "thermostats/thermostat.h"
#include "core/statistics_manager.h"
#include "strategies/observables/classical_spring_energy_strategy.h"

ClassicalObservable::ClassicalObservable(
    const std::shared_ptr<const VecArray>& coord,
    const std::shared_ptr<const VecArray>& prev_coord,
    const ThermostatContext& thermostat_ctx,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    int out_freq,
    const std::string& out_unit
) : Observable("classical", out_freq, out_unit),
    m_coord(coord),
    m_prev_coord(prev_coord),
    m_spring_energy_strategy(StatisticsManager::getInstance().createClassicalSpringEnergyStrategy()),
    m_thermostat_ctx(thermostat_ctx),
    m_spring_ctx(spring_ctx),
    m_box_ctx(box_ctx)
{
    m_is_nose_hoover = thermostat_ctx.thermostat_type.find("nose_hoover") != std::string::npos;
    //const bool is_nose_hoover = thermostat_ctx.thermostat_type == "nose_hoover" || thermostat_ctx.thermostat_type == "nose_hoover_np" || thermostat_ctx.thermostat_type == "nose_hoover_np_dim";

    initialize({ "cl_spring" });

    if (m_is_nose_hoover)
        initializeLabel("nh_energy");
}

ClassicalObservable::~ClassicalObservable() = default;

void ClassicalObservable::calculate() {
    calculateSpringEnergy();
    calculateThermostatEnergy();
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