#include "walkers/roulette_splitting_cycle_prob.h"

RouletteSplittingCycleProb::RouletteSplittingCycleProb(int nworlds, int local_rank, int walker_id, MPI_Comm& bead_world, std::mt19937& rand_gen) : 
    RouletteSplitting(nworlds, local_rank, walker_id, bead_world, rand_gen) {}

std::vector<double> RouletteSplittingCycleProb::evaluateImportance() {
    std::vector<double> importanceWeights(nworlds, 0.0);
    importanceWeights[0] = nworlds;
    return importanceWeights;
}