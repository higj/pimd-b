#include "thermostats/null_thermostat.h"

NullThermostat::NullThermostat(
    const ThermalContext& thermal_ctx,
    const NormalModesContext& nm_ctx,
    const std::shared_ptr<SystemState>& state
) : Thermostat(thermal_ctx, nm_ctx, state) {
    // Explicitly don't create a coupling
}