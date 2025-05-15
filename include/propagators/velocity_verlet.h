#pragma once

#include "propagators/propagator.h"

class Simulation;

class VelocityVerletPropagator : public Propagator {
public:
    VelocityVerletPropagator(Simulation& _sim);
    ~VelocityVerletPropagator() override = default;

    void step() override;
private:
    void momentStep();
    void coordsStep();
};