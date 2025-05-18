#pragma once

#include "common.h"

#include <memory>
#include <string>

class ExchangeState;
class ForceFieldManager;

/**
 * Parameters needed to calculate energy observables.
 */
struct EnergyObservableContext {
    std::shared_ptr<ExchangeState> exchange_state;
    std::shared_ptr<const dVec> coord;
    std::shared_ptr<const dVec> prev_coord;
    std::shared_ptr<ForceFieldManager> force_mgr;
    int natoms;
    int nbeads;
    int this_bead;
    double beta;
    double spring_constant;
    double box_size;
    bool bosonic;
    std::string ext_pot_name;
    std::string int_pot_name;
};