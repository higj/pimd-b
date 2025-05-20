#include "walkers/roulette_splitting_cycle_prob.h"
#include "bosonic_exchange.h"
#include <iostream>


RouletteSplittingCycleProb::RouletteSplittingCycleProb(int nworlds, int local_rank, int walker_id, MPI_Comm& bead_world, 
                                                       std::mt19937& rand_gen, BosonicExchangeBase& bosonic_exchange) : 
    RouletteSplitting(nworlds, local_rank, walker_id, bead_world, rand_gen), bosonic_exchange(bosonic_exchange) {}

double RouletteSplittingCycleProb::evaluateImportance() {
    double importance = bosonic_exchange.getCumulativeCycleProb();
    return importance;
}