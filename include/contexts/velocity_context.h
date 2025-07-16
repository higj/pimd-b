#pragma once

#include "common.h"
#include <memory>

/**
 * Data associated with the velocity of the beads.
 */
struct VelocityContext {
    std::shared_ptr<dVec> momenta;
    double mass;
};