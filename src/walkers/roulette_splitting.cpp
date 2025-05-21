#include "walkers/roulette_splitting.h"
#include <numeric>
#include <vector>
#include "params.h"
#include <algorithm>
#include <iostream>

RouletteSplitting::RouletteSplitting(Params& param_obj, int nworlds,int local_rank, int walker_id, MPI_Comm& bead_world, std::mt19937& rand_gen) : 
    nworlds(nworlds), local_rank(local_rank), walker_id(walker_id), 
    bead_world(bead_world), rand_gen(rand_gen), u_dist(0.0, 1.0) {
        getVariant(param_obj.sim["roullete_thershold"], roullete_thershold);
        getVariant(param_obj.sim["splitting_thershold"], splitting_thershold);
        getVariant(param_obj.sim["roullete_splitting_normalized_thershold"], roullete_splitting_normalized_thershold);
    }

void RouletteSplitting::communicate(dVec& coord, dVec& momenta) {
    
    MPI_Barrier(MPI_COMM_WORLD);

    int assigned_worlds[nworlds];
    for (int i = 0; i < nworlds; ++i) {
        assigned_worlds[i] = nworlds;
    }

    if (local_rank == 0) {        
        // Evaluate the importance weights
        importance_weight = evaluateImportance();
        std::vector<double> importance_weights(nworlds);
        MPI_Gather(&importance_weight, 1, MPI_DOUBLE, importance_weights.data(), 1, MPI_DOUBLE, 0, bead_world);
        std::vector<double> old_statistical_weights(nworlds);
        MPI_Gather(&statistical_weight, 1, MPI_DOUBLE, old_statistical_weights.data(), 1, MPI_DOUBLE, 0, bead_world);
        std::vector<double> new_statistical_weights(nworlds);
        MPI_Gather(&statistical_weight, 1, MPI_DOUBLE, new_statistical_weights.data(), 1, MPI_DOUBLE, 0, bead_world);

        if (walker_id == 0) {
            std::vector<int> ncopies(nworlds, 0);
            // Create a vector of ids for worlds participating, and a vector of weights
            double normalization = static_cast<double>(nworlds) / std::accumulate(importance_weights.begin(), importance_weights.end(), 0.0);
            std::vector<int> participating_worlds;
            std::vector<double> participating_weights;
            participating_weights.push_back(0.0);
            int n_participating_worlds = 0;
            for (int i = 0; i < nworlds; ++i) {
                bool allowed_to_die = ((importance_weights[i] * normalization < 1) && (importance_weights[i] < roullete_thershold));
                bool allowed_to_split = ((importance_weights[i] * normalization > 1) && (importance_weights[i] > splitting_thershold));
                if (allowed_to_die || allowed_to_split) {
                    participating_worlds.push_back(i);
                    participating_weights.push_back(importance_weights[i]);
                    n_participating_worlds++;
                } else {
                    ncopies[i] = 1;
                }
            }

            // Normalize the participating weights
            double sum = std::accumulate(participating_weights.begin(), participating_weights.end(), 0.0);
            for (double& v : participating_weights) {
                v /= sum;
            }
            
            // Create an array of random numbers
            std::vector<double> u_samples(n_participating_worlds);
            for (int i = 0; i < n_participating_worlds; ++i) {
                u_samples[i] = u_dist(rand_gen);
            }

            // Sample copies of the worlds
            double lower_limit = 0;
            double upper_limit = 0;
            for (int i = 0; i < n_participating_worlds; ++i) {
                lower_limit += participating_weights[i];
                upper_limit += participating_weights[i+1];
                for (int j = 0; j < n_participating_worlds; ++j) {
                    double u_sample = u_samples[j];
                    if ((u_sample > lower_limit) && (u_sample < upper_limit)) {
                        ncopies[participating_worlds[i]]++;
                    }
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
            for (int i = 0; i < n_participating_worlds; ++i) {
                new_statistical_weights[participating_worlds[i]] = 
                    old_statistical_weights[assigned_worlds[participating_worlds[i]]] / 
                    (n_participating_worlds * importance_weights[assigned_worlds[participating_worlds[i]]] / sum);
            }
        }
        // Broadcast the weights to all worlds
        MPI_Bcast(new_statistical_weights.data(), nworlds, MPI_DOUBLE, 0, bead_world);
        statistical_weight = new_statistical_weights[walker_id];
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