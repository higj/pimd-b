#pragma once

#include <vector>
#include <memory>

class Simulation;
class Coupling;

class Thermostat {
public:
    explicit Thermostat(Simulation& _sim, bool normal_modes);
    virtual ~Thermostat() = default;
    void step();
    virtual void momentaUpdate();
    virtual double getAdditionToH();
protected:
    Simulation& sim;   // Reference to the simulation object
    std::unique_ptr<Coupling> coupling;
};