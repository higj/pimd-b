#include "observables/bosonic_prob_all_observable.h"
#include "core/statistics_manager.h"
#include "strategies/observables/bosonic_probability_strategy.h"

BosonicProbAllObservable::BosonicProbAllObservable(
    int out_freq, const std::string& out_unit
) : Observable("prob_all", out_freq, out_unit),
m_strategy(StatisticsManager::getInstance().createBosonicProbabilityStrategy()) {
    initializeLabel("prob_all");
}

BosonicProbAllObservable::~BosonicProbAllObservable() = default;

void BosonicProbAllObservable::calculate() {
    quantities["prob_all"] = m_strategy->getLongestProbability();
}