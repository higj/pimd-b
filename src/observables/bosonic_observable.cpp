#include "observables/bosonic_observable.h"
/*
// Unquoting the following seems to speed up the simulation in case of MSVC compilation
namespace Units {
    const std::unordered_map<std::string, std::unordered_map<std::string, double>> UnitMap = {
        {"", {
            {"", 0}
        }}
    };
}
*/
#include "core/statistics_manager.h"
#include "strategies/observables/bosonic_probability_strategy.h"

BosonicObservable::BosonicObservable(int out_freq, const std::string& out_unit) : Observable("bosonic", out_freq, out_unit),
    m_bosonic_prob_strategy(StatisticsManager::getInstance().createBosonicProbabilityStrategy())
{
    initialize({"prob_dist", "prob_all", "sign" });
}

BosonicObservable::~BosonicObservable() = default;

void BosonicObservable::calculate()
{
    quantities["prob_dist"] = m_bosonic_prob_strategy->getDistinctProbability();
    quantities["prob_all"] = m_bosonic_prob_strategy->getLongestProbability();
    quantities["sign"] = m_bosonic_prob_strategy->getSign();
}
