#pragma once

#include "common.h"
#include "contexts/spring_context.h"
#include "contexts/box_context.h"

class ClassicalSpringEnergyStrategy {
public:
    virtual double calculateSpringEnergy(
        const dVec& coord,
        const dVec& prev_coord,
        const SpringContext& spring_ctx, 
        const BoxContext& box_ctx
    ) = 0;
    virtual ~ClassicalSpringEnergyStrategy() = default;
};