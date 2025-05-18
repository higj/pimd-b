#pragma once

#include <memory>
#include "common.h"

struct NormalModesContext {
    std::shared_ptr<const dVec> coord;
    std::shared_ptr<dVec> momenta;
    int natoms;
    int nbeads;
    int this_bead;
};