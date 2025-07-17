#pragma once

#include "interior_bosonic_probability_strategy.h"

#include <memory>

class BosonicExchangeBase;

class ExteriorBosonicProbabilityStrategy : public InteriorBosonicProbabilityStrategy {
public:
    ExteriorBosonicProbabilityStrategy(const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange);

    double getDistinctProbability() override;
    double getLongestProbability() override;

private:
    std::shared_ptr<BosonicExchangeBase> m_bosonic_exchange;
};
