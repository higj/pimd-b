#pragma once

class BosonicProbabilityStrategy {
public:
    virtual double getDistinctProbability() = 0;
    virtual double getLongestProbability() = 0;

    virtual ~BosonicProbabilityStrategy() = default;
};
