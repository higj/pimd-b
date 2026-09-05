#include "observables/bosonic_prob_dist_observable.h"
#include "core/statistics_manager.h"
#include "strategies/observables/bosonic_probability_strategy.h"

BosonicProbDistObservable::BosonicProbDistObservable(
    int out_freq, const std::string& out_unit
) : Observable("prob_dist", out_freq, out_unit),
m_strategy(StatisticsManager::getInstance().createBosonicProbabilityStrategy()) {
    initializeLabel("prob_dist");
}

BosonicProbDistObservable::~BosonicProbDistObservable() = default;

void BosonicProbDistObservable::calculate() {
    quantities["prob_dist"] = m_strategy->getDistinctProbability();
}