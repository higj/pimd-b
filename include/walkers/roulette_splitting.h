#pragma once

#include "walkers_communication_base.h"
#include <vector>
#include <random>
#include "mpi.h"

class RouletteSplitting : public WalkersCommunicationBase
{
public:
    explicit RouletteSplitting(int nworlds, int local_rank, int walker_id, MPI_Comm& bead_world, std::mt19937& rand_gen);
    ~RouletteSplitting() override = default;

    void communicate(dVec& coord, dVec& momenta) override;
    virtual std::vector<double> evaluateImportance() = 0;
protected:
    MPI_Comm& bead_world;
    int nworlds, local_rank, walker_id;
    std::mt19937& rand_gen;
    std::uniform_real_distribution<double> u_dist;
};