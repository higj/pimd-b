#pragma once

#include <memory>
#include <string>

class Thermostat;

struct ThermostatContext {
    std::shared_ptr<Thermostat> thermostat;
    std::string thermostat_type;
};