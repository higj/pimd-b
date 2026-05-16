#pragma once

#include "common.h"
#include "contexts/box_context.h"

namespace RingPolymerUtils
{
    double classicalSpringEnergy(const VecArray& coord, const VecArray& prev_coord, double spring_constant, const BoxContext& box_ctx);
}