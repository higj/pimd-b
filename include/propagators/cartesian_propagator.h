#pragma once

#include "propagators/propagator.h"

class CartesianPropagator : public Propagator {
public:
    CartesianPropagator(Params& param_obj, dVec& coord, dVec& momenta, dVec& forces);
    ~CartesianPropagator() override = default;

    void preForceStep() override;
    void postForceStep() override;
private:
    void momentStep();
    void coordsStep();
};