#pragma once

#include "roulette_splitting.h"
#include <vector>
#include "mpi.h"

class RouletteSplittingCycleProb : public RouletteSplitting
{
public:
    explicit RouletteSplittingCycleProb(int nworlds, int local_rank, int walker_id, MPI_Comm& bead_world, std::mt19937& rand_gen);
    ~RouletteSplittingCycleProb() override = default;

    std::vector<double> evaluateImportance() override;
};