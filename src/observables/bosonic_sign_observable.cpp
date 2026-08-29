#include "observables/bosonic_sign_observable.h"
#include "core/statistics_manager.h"
#include "strategies/observables/bosonic_probability_strategy.h"

BosonicSignObservable::BosonicSignObservable(
    int out_freq, const std::string& out_unit
) : Observable("sign", out_freq, out_unit),
m_strategy(StatisticsManager::getInstance().createBosonicProbabilityStrategy()) {
    initializeLabel("sign");
}

BosonicSignObservable::~BosonicSignObservable() = default;

void BosonicSignObservable::calculate() {
    quantities["sign"] = m_strategy->getSign();
}