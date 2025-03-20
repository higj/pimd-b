#include "RPC_simulation.h"
#include "mpi.h"
#include <numbers>
#include <math.h>

RPCSimulation::RPCSimulation(const int& rank, const int& nproc, Params& param_obj, unsigned int seed) : 
    Simulation(rank, nproc, param_obj, seed) {
    getVariant(param_obj.sim["nbeads_contructed"], nbeads_contructed);
    RPC_cutoff = std::get<double>(param_obj.interaction_pot["RPC_cutoff"]);
    beads_ratio = static_cast<double>(nbeads) / nbeads_contructed;

    // Allocate memory for transformation matrix T_{j',j}
    net_transformation_matrix = (double**)malloc(nbeads_contructed * sizeof(double*));
    for (int i = 0; i < nbeads_contructed; i++) {
        net_transformation_matrix[i] = (double*)malloc(nbeads * sizeof(double));
    }
    
    // Fill C_{j,k} according to Eq. 16-17 in Markland & Manolopoulos https://doi.org/10.1063/1.2953308
    double transformation_matrix[nbeads][nbeads];
    const double s_2 = sqrt(2);
    int l = ((nbeads % 2 == 0) ? (nbeads - 2) / 2 : (nbeads - 1) / 2);
    for (int j = 1; j < nbeads + 1; ++j) {
        transformation_matrix[j - 1][0] = 1;
        for (int k = 0; k < l; ++k) {
            transformation_matrix[j - 1][1 + k * 2] = s_2 * sin(-2 * j * (k+1) * M_PI / nbeads);
            transformation_matrix[j - 1][2 + k * 2] = s_2 * cos(2 * j * (k+1) * M_PI / nbeads);
        }
        if (nbeads % 2 == 0) {
            transformation_matrix[j - 1][nbeads - 1] = (j % 2 == 0 ? 1.0 : -1.0);
        }
    }

    // Fill T_{j',j} according to Eq. 22 in Markland & Manolopoulos
    for (int j_tag = 0; j_tag < nbeads_contructed; ++j_tag) {
        for (int j = 0; j < nbeads; ++j) {
            double new_val = 0;
            for (int k = 0; k < nbeads_contructed; ++k) {
                new_val += transformation_matrix[j_tag][k] * transformation_matrix[j][k];
            }
            net_transformation_matrix[j_tag][j] = new_val / nbeads;
        }
    }

    if (this_bead == 0) {
        MPI_Win_allocate_shared(nbeads*NDIM*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &ptcl_one_coords, &win_ptcl_one_coords);
        MPI_Win_allocate_shared(nbeads*NDIM*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &ptcl_two_coords, &win_ptcl_two_coords);
        MPI_Win_allocate_shared(nbeads_contructed*NDIM*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &constructed_forces, &win_constructed_forces);
    } else {
        int disp_unit;
        MPI_Aint size;
        MPI_Win_allocate_shared(0, sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &ptcl_one_coords, &win_ptcl_one_coords);
        MPI_Win_shared_query(win_ptcl_one_coords, 0, &size, &disp_unit, &ptcl_one_coords);
        MPI_Win_allocate_shared(0, sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &ptcl_two_coords, &win_ptcl_two_coords);
        MPI_Win_shared_query(win_ptcl_two_coords, 0, &size, &disp_unit, &ptcl_two_coords);
        MPI_Win_allocate_shared(0, sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &constructed_forces, &win_constructed_forces);
        MPI_Win_shared_query(win_constructed_forces, 0, &size, &disp_unit, &constructed_forces);
    }
}

RPCSimulation::~RPCSimulation() {
    for (int i = 0; i < nbeads_contructed; i++) {
        free(net_transformation_matrix[i]);
    }
    free(net_transformation_matrix);
}

void RPCSimulation::updatePhysicalForces(dVec& physical_force_arr) {
    // Calculate the external forces acting on the particles
    physical_force_arr = (-1.0) * ext_potential->gradV(coord);

    if (int_pot_cutoff != 0.0) {
        for (int ptcl_one = 0; ptcl_one < natoms; ++ptcl_one) {
            for (int ptcl_two = ptcl_one + 1; ptcl_two < natoms; ++ptcl_two) {
                // Get the vector distance between the two particles.
                // Here "diff" contains just one vector of dimension NDIM.
                dVec diff = getSeparation(ptcl_one, ptcl_two, MINIM);
                double norm = diff.norm();
                double summed_diff = 0;
                MPI_Allreduce(&norm, &summed_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                // If the average distance between the particles exceeds the cutoff length
                // then we assume the interaction is negligible and do not bother
                // calculating the force.
                // We use the convention that when cutoff < 0 then the interaction is
                // calculated for all distances.
                if (summed_diff / nbeads < int_pot_cutoff || int_pot_cutoff < 0.0) {
                    // If the average distance between particles is shorter than the RPC cutoff,
                    // calculate on all beads. Otherwise, only on contructed ring.
                    if (summed_diff / nbeads < RPC_cutoff) {
                        dVec force_on_one = (-1.0) * int_potential->gradV(diff);

                        for (int axis = 0; axis < NDIM; ++axis) {
                            physical_force_arr(ptcl_one, axis) += force_on_one(0, axis);
                            physical_force_arr(ptcl_two, axis) -= force_on_one(0, axis);
                        }
                    }
                    else {
                        updateForcesConstructedRing(physical_force_arr, ptcl_one, ptcl_two);
                    }
                }
            }
        }
    }
}

void RPCSimulation::updateForcesConstructedRing(dVec& physical_force_arr, int ptcl_one, int ptcl_two) {
    // Share positions of all beads of the two particles
    for (int axis = 0; axis < NDIM; ++axis) {
        ptcl_one_coords[this_bead * NDIM + axis] = coord(ptcl_one, axis);
        ptcl_two_coords[this_bead * NDIM + axis] = coord(ptcl_two, axis);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Calculate forces using contructed ring
    if (this_bead < nbeads_contructed) {
        dVec sep;
        for (int axis = 0; axis < NDIM; ++axis) {
            double coord_one = 0.;
            double coord_two = 0.;
            for (int j = 0; j < nbeads; ++j) {
                coord_one += net_transformation_matrix[this_bead][j] * ptcl_one_coords[j * NDIM + axis];
                coord_two += net_transformation_matrix[this_bead][j] * ptcl_two_coords[j * NDIM + axis];
            }
            double diff = coord_one - coord_two;
            if (pbc && MINIM)
                applyMinimumImage(diff, size);
            sep(0, axis) = diff;
        }

        dVec force_on_one = (-1.0) * int_potential->gradV(sep);
        for (int axis = 0; axis < NDIM; ++axis) {
            constructed_forces[this_bead * NDIM + axis] = force_on_one(0, axis);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Update the forces on the whole ring
    for (int j_tag = 0; j_tag < nbeads_contructed; ++j_tag) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double force_contribution = beads_ratio * net_transformation_matrix[j_tag][this_bead] * constructed_forces[j_tag * NDIM + axis];
            physical_force_arr(ptcl_one, axis) += force_contribution;
            physical_force_arr(ptcl_two, axis) -= force_contribution;
        }
    }
}
