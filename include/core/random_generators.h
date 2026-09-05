#pragma once

#include "random_mars.h"
#include <random>

class RandomGenerators {
public:
    explicit RandomGenerators(unsigned int seed = std::random_device{}());

    double uniform();                        // [0, 1)
    double uniform(double a, double b);      // [a, b)
    double normal(double mean = 0.0, double stddev = 1.0);
    double gaussian();
    int integer(int min, int max);

    std::mt19937& engine();

private:
    std::mt19937 m_rng;
    RanMars m_mars_gen;
    std::uniform_real_distribution<double> m_unit_dist;
};
