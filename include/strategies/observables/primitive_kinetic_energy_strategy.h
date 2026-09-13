#pragma once

#include "common.h"
#include "contexts/spring_context.h"
#include "contexts/box_context.h"
#include "contexts/bead_context.h"

class PrimitiveKineticEnergyStrategy {
public:
    virtual double calculateSpringContribution(
        const VecArray& coord,
        const VecArray& prev_coord,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        const BeadContext& bead_ctx
    ) = 0;
    virtual ~PrimitiveKineticEnergyStrategy() = default;
};