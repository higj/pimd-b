#pragma once

#include "contexts/spring_context.h"

#include <memory>

class SystemState;

class NormalModesMomentaStrategy {
public:
    virtual void momentaExternalForces(
        const std::shared_ptr<SystemState>& state,
        const SpringContext& spring_ctx,
        double dt
    ) = 0;
    virtual ~NormalModesMomentaStrategy() = default;
};