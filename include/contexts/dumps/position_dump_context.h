#pragma once

#include "common.h"
#include <memory>

/**
 * Parameters required for the positions dump.
 */
struct PositionDumpContext {
    std::shared_ptr<const dVec> coord;
    int natoms;
    int this_bead;
};