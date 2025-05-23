#pragma once

#include "common.h"

#include <memory>
#include <string>

class ExchangeState;
class Thermostat;

/**
 * Parameters needed to calculate classical observables.
 */
struct ClassicalObservableContext {
    std::shared_ptr<const dVec> coord;
    std::shared_ptr<const dVec> prev_coord;
    std::shared_ptr<const dVec> momenta;
    std::shared_ptr<ExchangeState> exchange_state;
    std::shared_ptr<Thermostat> thermostat;
    int natoms;
    int nbeads;
    int this_bead;
    double beta;
    double mass;
    double spring_constant;
    double box_size;
    bool pbc;
    bool bosonic;
    std::string thermostat_type;
};