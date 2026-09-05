#pragma once

#include <memory>
#include <string>

class Thermostat;

struct ThermostatContext {
    std::shared_ptr<Thermostat> thermostat;
    std::string thermostat_type;

    [[nodiscard]] bool isEnabled() const { return !thermostat_type.empty() && thermostat_type != "none"; }
};