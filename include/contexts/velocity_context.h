#pragma once

#include "common.h"
#include <memory>

/**
 * Data associated with the velocity of the beads.
 */
struct VelocityContext {
    std::shared_ptr<VecArray> momenta;
    double mass;
};