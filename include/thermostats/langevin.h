#pragma once

#include <vector>
#include <memory>
#include "thermostats/thermostat.h"

class Simulation;
class Coupling;

class LangevinThermostat : public Thermostat {
public:
    LangevinThermostat(Simulation& _sim, bool normal_modes);
    ~LangevinThermostat() override = default;

    void momentaUpdate() override;

protected:
    double friction_coefficient, noise_coefficient;
};