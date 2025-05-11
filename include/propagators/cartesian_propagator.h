#pragma once

#include "propagators/propagator.h"

class CartesianPropagator : public Propagator {
public:
    CartesianPropagator(Params& param_obj, dVec& coord, dVec& momenta, dVec& forces);
    ~CartesianPropagator() override = default;

    virtual void preForceStep();
    virtual void postForceStep();
private:
    void momentStep();
    void coordsStep();
};