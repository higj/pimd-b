#pragma once

#include "common.h"

class Simulation;
class Params;

class Propagator {
public:
    explicit Propagator(Simulation& _sim, Params& param_obj, dVec& coord, dVec& momenta, dVec& forces);
    virtual ~Propagator() = default;
    
    virtual void preForceStep() = 0;
    virtual void postForceStep() = 0;

protected:
    Simulation& sim; // Reference to the simulation object
    dVec& coord;
    dVec& momenta;
    dVec& forces;
    int natoms;
    double mass;
    double dt;
};