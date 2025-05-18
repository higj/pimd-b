#pragma once

//#include <memory>
#include "thermostats/thermostat.h"
#include "contexts/thermostats/langevin_thermostat_context.h"

class ThermostatContext;

class LangevinThermostat final : public Thermostat {
public:
    LangevinThermostat(const ThermostatContext& context, const LangevinThermostatContext& langevin_context);
    ~LangevinThermostat() override = default;

    void momentaUpdate() override;
private:
    double m_friction_coefficient;
    double m_noise_coefficient;
    LangevinThermostatContext m_langevin_context;
};