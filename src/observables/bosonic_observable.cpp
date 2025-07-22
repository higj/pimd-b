#include "observables/bosonic_observable.h"
#include "core/statistics_manager.h"
#include "strategies/observables/bosonic_probability_strategy.h"

BosonicObservable::BosonicObservable(int out_freq, const std::string& out_unit) : Observable(out_freq, out_unit),
    m_bosonic_prob_strategy(StatisticsManager::getInstance().createBosonicProbabilityStrategy())
{
    initialize({"prob_dist", "prob_all"});
}

BosonicObservable::~BosonicObservable() = default;

void BosonicObservable::calculate()
{
    quantities["prob_dist"] = m_bosonic_prob_strategy->getDistinctProbability();
    quantities["prob_all"] = m_bosonic_prob_strategy->getLongestProbability();
}
