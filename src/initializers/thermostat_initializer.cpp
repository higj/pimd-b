#include "initializers/thermostat_initializer.h"
#include "core/simulation_config.h"
#include "thermostats.h"

ThermostatInitializer::ThermostatInitializer(
    const std::shared_ptr<SimulationConfig>& config,
    const ThermalContext& thermal_ctx,
    const NormalModesContext& nm_ctx
) :
    m_mass(config->mass),
    m_dt(config->dt),
    m_thermostat_type(config->thermostat_type),
    m_thermostat_params(config->thermostat_params),
    m_thermal_ctx(thermal_ctx),
    m_nm_ctx(nm_ctx)
{
}

std::shared_ptr<Thermostat> ThermostatInitializer::createThermostat(
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<RandomGenerators>& rng
) const {
    if (m_thermostat_type == "langevin") {
        double gamma = std::get<double>(m_thermostat_params.at("gamma"));

        return std::make_shared<LangevinThermostat>(
            m_thermal_ctx,
            m_nm_ctx,
            state,
            rng,
            gamma,
            m_dt,
            m_mass
        );
    }

    int nchains = std::get<int>(m_thermostat_params.at("nchains"));

    if (m_thermostat_type == "nose_hoover") {
        return std::make_shared<NoseHooverThermostat>(m_thermal_ctx, m_nm_ctx, state, nchains, m_dt, m_mass);
    }

    if (m_thermostat_type == "nose_hoover_np") {
        return std::make_shared<NoseHooverNpThermostat>(m_thermal_ctx, m_nm_ctx, state, nchains, m_dt, m_mass);
    }

    if (m_thermostat_type == "nose_hoover_np_dim") {
        return std::make_shared<NoseHooverNpDimThermostat>(m_thermal_ctx, m_nm_ctx, state, nchains, m_dt, m_mass);
    }

    if (m_thermostat_type == "none") {
        return std::make_shared<Thermostat>(m_thermal_ctx, m_nm_ctx, state);
    }

    /// TODO: Raise an exception if the thermostat type is unknown?
    return {};
}