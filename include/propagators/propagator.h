#pragma once

#include "common.h"

class Params;

class Propagator {
public:
    explicit Propagator(Params& param_obj, dVec& coord, dVec& momenta, dVec& forces);
    virtual ~Propagator() = default;
    
    virtual void preForceStep() = 0;
    virtual void postForceStep() = 0;

protected:
    dVec& coord;
    dVec& momenta;
    dVec& forces;
    int natoms;
    double mass;
    double dt;
};