#include "observables/bosonic.h"
#include "simulation.h"
#include <ranges>
#include "mpi.h"

/**
 * @brief Bosonic observable class constructor.
 */
BosonicObservable::BosonicObservable(const Simulation& _sim, int _freq, const std::string& _out_unit) :
    Observable(_sim, _freq, _out_unit) {
    initialize({ "prob_dist", "prob_all" });
}

/**
 * @brief Calculates quantities pertaining to bosonic exchange.
 */
void BosonicObservable::calculate() {
    if (sim.this_bead == 0 && sim.bosonic) {
        quantities["prob_dist"] = sim.bosonic_exchange->getDistinctProbability();
        quantities["prob_all"] = sim.bosonic_exchange->getLongestProbability();
    }
}