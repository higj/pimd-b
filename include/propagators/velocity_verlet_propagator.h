#pragma once

#include "propagators/propagator.h"

class VelocityVerletPropagator final : public Propagator {
public:
    VelocityVerletPropagator(const PropagatorContext& context);
    ~VelocityVerletPropagator() override = default;

    void step() override;
};