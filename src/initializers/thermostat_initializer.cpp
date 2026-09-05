#include "initializers/thermostat_initializer.h"
#include "thermostats.h"

ThermostatInitializer::ThermostatInitializer(
    const std::shared_ptr<SimulationConfig>& config,
    const ThermalContext& thermal_ctx,
    const NormalModesContext& nm_ctx
) :
    m_mass(config->mass),
    m_dt(config->dt),
    m_thermostat_config(config->thermostat),
    m_thermal_ctx(thermal_ctx),
    m_nm_ctx(nm_ctx)
{
}

std::shared_ptr<Thermostat> ThermostatInitializer::createThermostat(
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<RandomGenerators>& rng
) const {
    if (m_thermostat_config.type == "none") {
        return std::make_shared<NullThermostat>(m_thermal_ctx, m_nm_ctx, state);
    }

    if (m_thermostat_config.type == "langevin") {
        return std::make_shared<LangevinThermostat>(
            m_thermal_ctx,
            m_nm_ctx,
            state,
            rng,
            m_thermostat_config.gamma,
            m_dt,
            m_mass
        );
    }

    if (m_thermostat_config.is_nose_hoover) {
        if (m_thermostat_config.type == "nose_hoover") {
            return std::make_shared<NoseHooverThermostat>(m_thermal_ctx, m_nm_ctx, state, m_thermostat_config.nchains, m_dt, m_mass);
        }

        if (m_thermostat_config.type == "nose_hoover_np") {
            return std::make_shared<NoseHooverNpThermostat>(m_thermal_ctx, m_nm_ctx, state, m_thermostat_config.nchains, m_dt, m_mass);
        }

        if (m_thermostat_config.type == "nose_hoover_np_dim") {
            return std::make_shared<NoseHooverNpDimThermostat>(m_thermal_ctx, m_nm_ctx, state, m_thermostat_config.nchains, m_dt, m_mass);
        }
    }

    /// TODO: Raise an exception if the thermostat type is unknown?
    return {};
}