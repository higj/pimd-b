#pragma once

#include "propagators/propagator.h"

class Simulation;

class VelocityVerletPropagator : public Propagator {
public:
    VelocityVerletPropagator(Simulation& _sim, Params& param_obj, dVec& coord, dVec& momenta, dVec& forces);
    ~VelocityVerletPropagator() override = default;

    void preForceStep() override;
    void postForceStep() override;
private:
    void momentStep();
    void coordsStep();
};