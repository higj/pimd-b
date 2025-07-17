#include "observables/bosonic_observable.h"
#include "bosonic_exchange/bosonic_exchange_base.h"
#include "core/statistics_manager.h"

BosonicObservable::BosonicObservable(
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
    const BeadContext& bead_ctx,
    bool bosonic,
    int out_freq,
    const std::string& out_unit
) : Observable(out_freq, out_unit),
    m_bosonic_prob_strategy(StatisticsManager::createBosonicProbabilityStrategy(bosonic_exchange, bead_ctx, bosonic))
{
    initialize({"prob_dist", "prob_all"});
}

void BosonicObservable::calculate()
{
    quantities["prob_dist"] = m_bosonic_prob_strategy->getDistinctProbability();
    quantities["prob_all"] = m_bosonic_prob_strategy->getLongestProbability();
}
