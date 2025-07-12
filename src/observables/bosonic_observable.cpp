#include "observables/bosonic_observable.h"
#include "core/exchange_state.h"
#include "bosonic_exchange/bosonic_exchange_base.h"

BosonicObservable::BosonicObservable(
    const std::shared_ptr<ExchangeState>& exchange_state,
    int out_freq,
    const std::string& out_unit
) : Observable(out_freq, out_unit),
    m_exchange_state(exchange_state)
{
    initialize({"prob_dist", "prob_all"});
}

void BosonicObservable::calculate()
{
    if (m_exchange_state->is_first_bead)
    {
        quantities["prob_dist"] = m_exchange_state->bosonic_exchange->getDistinctProbability();
        quantities["prob_all"] = m_exchange_state->bosonic_exchange->getLongestProbability();
    }
}
