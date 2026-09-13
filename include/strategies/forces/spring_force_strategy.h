#pragma once

#include "core/system_state.h"
#include "contexts/spring_context.h"
#include "contexts/box_context.h"

class SpringForceStrategy {
public:
    virtual void updateSpringForces(SystemState& state, const SpringContext& spring_ctx, const BoxContext& box_ctx) = 0;
    virtual ~SpringForceStrategy() = default;
};