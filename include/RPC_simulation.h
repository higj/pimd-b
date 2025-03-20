#pragma once

#include "simulation.h"
#include "mpi.h"

class RPCSimulation : public Simulation
{
public:
    RPCSimulation(const int& rank, const int& nproc, Params& param_obj, unsigned int seed = static_cast<unsigned int>(time(nullptr)));
    ~RPCSimulation();
    void updatePhysicalForces(dVec& physical_force_arr) override;
private:
    void updateForcesConstructedRing(dVec& physical_force_arr, int ptcl_one, int ptcl_two);
    double RPC_cutoff;
    int nbeads_contructed;
    double **net_transformation_matrix;
    double *ptcl_one_coords, *ptcl_two_coords, *constructed_forces;
    MPI_Win win_ptcl_one_coords, win_ptcl_two_coords, win_constructed_forces;     // Window objects associated with the arrays
    double beads_ratio;
};
