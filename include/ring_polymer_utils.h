#pragma once

#include "common.h"

namespace RingPolymerUtils
{
    double classicalSpringEnergy(const dVec& coord, const dVec& prev_coord, double spring_constant, bool minimum_image, double box_size);
}