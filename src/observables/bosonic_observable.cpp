#include "observables/bosonic_observable.h"
#include "core/exchange_state.h"
#include "bosonic_exchange/bosonic_exchange_base.h"

/**
 * @brief Bosonic observable class constructor.
 */
BosonicObservable::BosonicObservable(const BosonicObservableContext& obs_context, int out_freq, const std::string& out_unit) :
    Observable(out_freq, out_unit), m_context(obs_context)
{
    initialize({ "prob_dist", "prob_all" });
}

/**
 * @brief Calculates quantities pertaining to bosonic exchange.
 */
void BosonicObservable::calculate() {
    if (m_context.this_bead == 0) {
        quantities["prob_dist"] = m_context.exchange_state->bosonic_exchange->getDistinctProbability();
        quantities["prob_all"] = m_context.exchange_state->bosonic_exchange->getLongestProbability();
    }
}