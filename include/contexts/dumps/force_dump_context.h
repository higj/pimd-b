#pragma once

#include "common.h"
#include <memory>

class SystemState;
/**
 * Parameters required for the forces dump.
 */
struct ForceDumpContext {
    std::shared_ptr<SystemState> state;
};