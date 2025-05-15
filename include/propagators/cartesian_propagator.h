#pragma once

#include "propagators/propagator.h"

class VelocityVerletPropagator : public Propagator {
public:
    VelocityVerletPropagator(Params& param_obj, dVec& coord, dVec& momenta, dVec& forces);
    ~VelocityVerletPropagator() override = default;

    void preForceStep() override;
    void postForceStep() override;
private:
    void momentStep();
    void coordsStep();
};