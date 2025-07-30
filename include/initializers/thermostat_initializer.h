#pragma once

#include "common.h"
#include "contexts/thermal_context.h"
#include "contexts/normal_modes_context.h"

#include <memory>
#include <string>

struct SimulationConfig;
class Thermostat;
class SystemState;
class RandomGenerators;

class ThermostatInitializer {
public:
    ThermostatInitializer(
        const std::shared_ptr<SimulationConfig>& config,
        const ThermalContext& thermal_ctx,
        const NormalModesContext& nm_ctx
    );

    [[nodiscard]] std::shared_ptr<Thermostat> createThermostat(
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<RandomGenerators>& rng
    ) const;

private:
    double m_mass;
    double m_dt;
    std::string m_thermostat_type;
    VariantMap m_thermostat_params;
    ThermalContext m_thermal_ctx;
    NormalModesContext m_nm_ctx;
};
