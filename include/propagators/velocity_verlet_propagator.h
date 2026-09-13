#pragma once

#include "propagators/propagator.h"

class VelocityVerletPropagator final : public Propagator
{
public:
    VelocityVerletPropagator(
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<ForceManager>& force_mgr,
        const BoxContext& box_ctx,
        const SpringContext& spring_ctx,
        double mass,
        double dt
    );
    ~VelocityVerletPropagator() override = default;

    void step() override;

    void momentStep() const;
    void coordsStep() const;
};
