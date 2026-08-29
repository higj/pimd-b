#pragma once

#include <cstddef>

enum class CacheKey : std::size_t {
    CL_KINETIC_RAW = 0,
    // add further keys here as the refactor proceeds
    COUNT  // sentinel — always last
};