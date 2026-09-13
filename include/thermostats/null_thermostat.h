#pragma once

#include "thermostats/thermostat.h"

// Represents the NVE (no thermostat) case
class NullThermostat : public Thermostat {
public:
    explicit NullThermostat(
        const ThermalContext& thermal_ctx,
        const NormalModesContext& nm_ctx,
        const std::shared_ptr<SystemState>& state
    );

    // Override to provide a no-op step that doesn't create coupling objects
    void step() override {} // Literally do nothing

    void momentaUpdate() override {} // Explicitly empty
    double getAdditionToH() override { return 0.0; }
};