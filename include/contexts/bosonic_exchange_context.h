#pragma once

#include "common.h"
#include <memory>

struct BosonicExchangeContext {
    int nbosons;
    int nbeads;
    double spring_constant;
    double beta_half_k;
    double beta;
    double thermo_beta;
    std::shared_ptr<const dVec> x;
    std::shared_ptr<const dVec> x_prev;
    std::shared_ptr<const dVec> x_next;
    bool pbc;
    double box_size;
    int this_bead;
};