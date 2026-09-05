#pragma once

#include <cstddef>

enum class CacheKey : std::size_t {
    CL_KINETIC_RAW = 0,      // classical kinetic energy, pre /nbeads, internal units
    EXT_POT_RAW = 1,         // ext_potential->V(coord), pre /nbeads, internal units
    INT_POT_RAW = 2,         // pair-loop sum of int_potential->V(diff), pre /nbeads
    VIRIAL_RAW = 3,          // -sum_i(r_i * F_i), pre *0.5/nbeads
    SC_TOTAL_POTENTIAL = 4,  // ext + int pair potential (raw)
    SC_FORCE_SQUARED = 5,    // sum |gradV|^2 (raw)
    SC_VIRIAL = 6,           // GSF virial sum (raw)
    COUNT = 7                // sentinel - always last
};