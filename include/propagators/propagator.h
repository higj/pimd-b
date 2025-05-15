#pragma once

#include "common.h"

class Simulation;

class Propagator {
public:
    explicit Propagator(Simulation& _sim);
    virtual ~Propagator() = default;
    
    virtual void preForceStep() = 0;
    virtual void postForceStep() = 0;

protected:
    Simulation& sim; // Reference to the simulation object
};