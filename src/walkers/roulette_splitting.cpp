#include "walkers/roulette_splitting.h"
#include <numeric>
#include <vector>
#include <iostream>

RouletteSplitting::RouletteSplitting(int nworlds,int local_rank, int walker_id, MPI_Comm& bead_world, std::mt19937& rand_gen) : 
    nworlds(nworlds), local_rank(local_rank), walker_id(walker_id), bead_world(bead_world), rand_gen(rand_gen), u_dist(0.0, 1.0) {}

void RouletteSplitting::communicate(dVec& coord, dVec& momenta) {
    
    MPI_Barrier(MPI_COMM_WORLD);

    int assigned_worlds[nworlds];
    for (int i = 0; i < nworlds; ++i) {
        assigned_worlds[i] = nworlds;
    }

    if (local_rank == 0) {
        int ncopies[nworlds];
        // Evaluate the importance weights
        importance_weight = evaluateImportance();
        std::vector<double> impotance_weights(nworlds);
        MPI_Gather(&importance_weight, 1, MPI_DOUBLE, impotance_weights.data(), 1, MPI_DOUBLE, 0, bead_world);
        std::vector<double> statistical_weights(nworlds);
        MPI_Gather(&statistical_weight, 1, MPI_DOUBLE, statistical_weights.data(), 1, MPI_DOUBLE, 0, bead_world);

        if (walker_id == 0) {
            // Normalize the importance weights
            double sum = std::accumulate(impotance_weights.begin(), impotance_weights.end(), 0.0);
            if (sum != 0.0) {
                double scale = static_cast<double>(nworlds) / sum;
                for (double& v : impotance_weights) {
                    v *= scale;
                }
            }

            int current_nworlds = 0;
            // Assign the number of copies of each world such that the total number of copies is equal to nworlds
            while (current_nworlds != nworlds) {
                current_nworlds = 0;
                for (int i = 0; i < nworlds; ++i) {
                    double u_sample = u_dist(rand_gen);
                    if (impotance_weights[i] < 1) {
                        ncopies[i] = (u_sample < impotance_weights[i]) ? 1 : 0;
                    }
                    else {
                        int floor_val = static_cast<int>(impotance_weights[i]);
                        ncopies[i] = (u_sample < impotance_weights[i] - floor_val) ? floor_val + 1 : floor_val;
                    }
                    current_nworlds += ncopies[i];
                }
            }
            // Assign a copy to each world that has at least one copy
            for (int i = 0; i < nworlds; ++i) {
                if (ncopies[i] > 0) {
                    assigned_worlds[i] = i;
                }
            }
            // Assign remaining copies
            for (int i = 0; i < nworlds; ++i) {
                if (ncopies[i] > 1) {
                    for (int n = 0; n < ncopies[i] -1; ++n) {
                        for (int j = 0; j < nworlds; ++j) {
                            if (assigned_worlds[j] == nworlds) {
                                assigned_worlds[j] = i;
                                break;
                            }
                        }
                    }
                }           
            }
            // Adjust statistical weights
            for (int i = 0; i < nworlds; ++i) {
                if ((impotance_weights[i] < 1) && (assigned_worlds[i] == i)) {
                    statistical_weights[i] *= 1 / impotance_weights[i];
                } else {
                    statistical_weights[i] = statistical_weights[assigned_worlds[i]] / ncopies[assigned_worlds[i]];
                }
            }
        }
        // Broadcast the statistical weights to all worlds
        MPI_Bcast(statistical_weights.data(), nworlds, MPI_DOUBLE, 0, bead_world);        
    }
    // Broadcast the assigned worlds to all worlds
    MPI_Bcast(assigned_worlds, nworlds, MPI_INT, 0, MPI_COMM_WORLD);

    // Broadcast the assigned copies to all worlds
    const int size = coord.size();
    for (int i = 0; i < nworlds; ++i) {
        if (assigned_worlds[i] != i) {
            // If in the world that is being copy, send coords and momenta to world i
            if (walker_id == assigned_worlds[i]) {
                MPI_Send(coord.data(), size, MPI_DOUBLE, i, 0, bead_world);
                MPI_Send(momenta.data(), size, MPI_DOUBLE, i, 1, bead_world);

            }
            // If in world i, receive coords and momenta from the world that is being copied
            else if (walker_id == i) {
                dVec new_coord(size);
                dVec new_momenta(size);
                MPI_Recv(new_coord.data(), size, MPI_DOUBLE, assigned_worlds[i], 0, bead_world, MPI_STATUS_IGNORE);
                MPI_Recv(new_momenta.data(), size, MPI_DOUBLE, assigned_worlds[i], 1, bead_world, MPI_STATUS_IGNORE);
                for (int ptcl_idx = 0; ptcl_idx < size/NDIM; ++ptcl_idx) {
                    for (int axis = 0; axis < NDIM; ++axis) {
                        coord(ptcl_idx, axis) = new_coord(ptcl_idx, axis);
                        momenta(ptcl_idx, axis) = new_momenta(ptcl_idx, axis);
                    }
                }
            }
        }
    }
}