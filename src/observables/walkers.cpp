#include "observables/walkers.h"
#include "walkers_communication.h"
#include <ranges>

/**
 * @brief Bosonic observable class constructor.
 */
WalkersObservable::WalkersObservable(Params& param_obj, int _freq, const std::string& _out_unit, 
                                     int this_bead, WalkersCommunicationBase& walker_communication) :
    Observable(param_obj, _freq, _out_unit, this_bead), walker_communication(walker_communication) {
    initialize({ "importance_weight", "statistical_weight" });
}

/**
 * @brief Calculates quantities pertaining to bosonic exchange.
 */
void WalkersObservable::calculate() {
    if (this_bead == 0) {
        quantities["importance_weight"] = walker_communication.importance_weight;
        quantities["statistical_weight"] = walker_communication.statistical_weight;
    }
}