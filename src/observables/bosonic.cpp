#include "observables/bosonic.h"
#include "bosonic_exchange.h"
#include <ranges>

/**
 * @brief Bosonic observable class constructor.
 */
BosonicObservable::BosonicObservable(Params& param_obj, int _freq, const std::string& _out_unit, 
                                     int this_bead, BosonicExchangeBase& bosonic_exchange) :
    Observable(param_obj, _freq, _out_unit, this_bead), bosonic_exchange(bosonic_exchange) {
    getVariant(param_obj.sim["bosonic"], bosonic);
    initialize({ "prob_dist", "prob_all" });
}

/**
 * @brief Calculates quantities pertaining to bosonic exchange.
 */
void BosonicObservable::calculate() {
    if (this_bead == 0 && bosonic) {
        quantities["prob_dist"] = bosonic_exchange.getDistinctProbability();
        quantities["prob_all"] = bosonic_exchange.getLongestProbability();
    }
}