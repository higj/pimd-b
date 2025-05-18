#pragma once

#include "common.h"
#include <memory>

/**
 * Parameters required for the velocities dump.
 */
struct VelocityDumpContext {
    std::shared_ptr<const dVec> momenta;
    int natoms;
    int this_bead;
    double mass;
};