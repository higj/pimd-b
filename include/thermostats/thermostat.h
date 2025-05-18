#pragma once

#include "contexts/thermostats/thermostat_context.h"

#include <memory>

class Coupling;

class Thermostat {
public:
    explicit Thermostat(const ThermostatContext& context);
    virtual ~Thermostat() = default;

    void step();
    virtual void momentaUpdate();
    virtual double getAdditionToH();
protected:
    ThermostatContext m_context;
    std::unique_ptr<Coupling> coupling;
};