#pragma once

#include "common.h"
#include "contexts/box_context.h"

namespace RingPolymerUtils
{
    double classicalSpringEnergy(const dVec& coord, const dVec& prev_coord, double spring_constant, const BoxContext& box_ctx);
}