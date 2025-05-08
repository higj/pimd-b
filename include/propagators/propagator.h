#pragma once

#include "common.h"

class Simulation;

class Propagator {
public:
    explicit Propagator(Simulation& _sim);
    virtual ~Propagator() = default;
    
    virtual void step() = 0;
    void momentStep();
    void coordsStep();

protected:
    Simulation& sim; // Reference to the simulation object
};