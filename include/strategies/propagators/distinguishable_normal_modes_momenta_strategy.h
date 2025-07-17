#pragma once

#include "normal_modes_momenta_strategy.h"

class DistinguishableNormalModesMomentaStrategy : public NormalModesMomentaStrategy {
public:
    void momentaExternalForces(
        const std::shared_ptr<SystemState>& state,
        const SpringContext& spring_ctx,
        double dt
    ) override;
};
