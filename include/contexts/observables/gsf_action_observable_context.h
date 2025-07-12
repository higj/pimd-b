#pragma once

#include "common.h"

class ForceManager;

/**
 * Parameters needed for the GSF action observables.
 */
struct GSFActionObservableContext {
    std::shared_ptr<const dVec> coord;
    std::shared_ptr<ForceManager> force_mgr;
    int natoms;
    int nbeads;
    int this_bead;
    double beta;
    double spring_constant;
};