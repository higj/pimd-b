#pragma once

#include <cstddef>

enum class CacheKey : std::size_t {
    CL_KINETIC_RAW = 0,  // existing
    EXT_POT_RAW = 1,     // ext_potential->V(coord), pre /nbeads, internal units
    INT_POT_RAW = 2,     // pair-loop sum of int_potential->V(diff), pre /nbeads
    VIRIAL_RAW = 3,      // -sum_i(r_i * F_i), pre *0.5/nbeads
    COUNT = 4            // sentinel - always last
};