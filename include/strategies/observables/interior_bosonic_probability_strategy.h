#pragma once

#include "bosonic_probability_strategy.h"

class InteriorBosonicProbabilityStrategy : public BosonicProbabilityStrategy {
public:
    double getDistinctProbability() override { return 0.0; }
    double getLongestProbability() override { return 0.0; }
};
