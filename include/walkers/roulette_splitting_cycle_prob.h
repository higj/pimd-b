#pragma once

#include "roulette_splitting.h"
#include <vector>
#include "mpi.h"

class BosonicExchangeBase;

class RouletteSplittingCycleProb : public RouletteSplitting
{
public:
    explicit RouletteSplittingCycleProb(int nworlds, int local_rank, int walker_id, MPI_Comm& bead_world, 
                                        std::mt19937& rand_gen, BosonicExchangeBase& bosonic_exchange);
    ~RouletteSplittingCycleProb() override = default;

    double evaluateImportance() override;
private:
    BosonicExchangeBase& bosonic_exchange;
};