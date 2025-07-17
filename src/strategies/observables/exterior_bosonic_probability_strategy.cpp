#include "strategies/observables/exterior_bosonic_probability_strategy.h"
#include "bosonic_exchange/bosonic_exchange_base.h"

ExteriorBosonicProbabilityStrategy::ExteriorBosonicProbabilityStrategy(
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange
) : m_bosonic_exchange(bosonic_exchange)
{
}

double ExteriorBosonicProbabilityStrategy::getDistinctProbability() {
    return m_bosonic_exchange->getDistinctProbability();
}


double ExteriorBosonicProbabilityStrategy::getLongestProbability() {
    return m_bosonic_exchange->getLongestProbability();
}