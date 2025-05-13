#pragma once

#include <vector>
#include <memory>
#include "thermostats/thermostat.h"
#include "random_mars.h"

class Simulation;
class Coupling;

class LangevinThermostat : public Thermostat {
public:
    LangevinThermostat(Coupling& coupling, Params& param_obj, int rank);
    ~LangevinThermostat() override = default;

    void momentaUpdate() override;

private:
    double friction_coefficient, noise_coefficient;
    std::unique_ptr<RanMars> mars_gen;
};