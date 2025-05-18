#pragma once

#include "common.h"
#include <memory>

class RandomGenerators;

/**
 * Parameters needed for the Box class.
 */
struct BoxContext {
    std::shared_ptr<dVec> coord;
    //std::shared_ptr<dVec> momenta;
    std::shared_ptr<RandomGenerators> rng;
    int natoms;
    double box_size;
    //bool mic_spring;
    std::string init_pos_type;
    //std::string init_vel_type;
};