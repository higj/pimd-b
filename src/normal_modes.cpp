#include <numbers>
#include "common.h"
#include "normal_modes.h"
#include "simulation.h"

NormalModes::NormalModes(Simulation& _sim) :
    sim(_sim),
    axis_stride(_sim.natoms*_sim.nbeads),
    atom_stride(_sim.nbeads),
    cart_to_nm_mat_row(_sim.nbeads),
    nm_to_cart_mat_row(_sim.nbeads)
{
    // Allocate shared memory
    if (sim.this_bead == 0) {
        MPI_Win_allocate_shared(sim.natoms*sim.nbeads*NDIM*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_coord_cartesian, &win_coord_cartesian);
        MPI_Win_allocate_shared(sim.natoms*sim.nbeads*NDIM*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_coord_nm, &win_coord_nm);
        MPI_Win_allocate_shared(sim.natoms*sim.nbeads*NDIM*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_momenta_cartesian, &win_momenta_cartesian);
        MPI_Win_allocate_shared(sim.natoms*sim.nbeads*NDIM*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_momenta_nm, &win_momenta_nm);
    } else {
        int disp_unit;
        MPI_Aint size;
        MPI_Win_allocate_shared(0, sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_coord_cartesian, &win_coord_cartesian);
        MPI_Win_allocate_shared(0, sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_coord_nm, &win_coord_nm);
        MPI_Win_allocate_shared(0, sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_momenta_cartesian, &win_momenta_cartesian);
        MPI_Win_allocate_shared(0, sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_momenta_nm, &win_momenta_nm);
        MPI_Win_shared_query(win_coord_cartesian, 0, &size, &disp_unit, &arr_coord_cartesian);
        MPI_Win_shared_query(win_coord_nm, 0, &size, &disp_unit, &arr_coord_nm);
        MPI_Win_shared_query(win_momenta_cartesian, 0, &size, &disp_unit, &arr_momenta_cartesian);
        MPI_Win_shared_query(win_momenta_nm, 0, &size, &disp_unit, &arr_momenta_nm);
    }
    
    // Frequencies
    double freq = 2 * sim.omega_p * sin(sim.this_bead * std::numbers::pi / sim.nbeads);
    
    // Cartesian-to-nm transformation matrix (one row because parallelized)
    double pref;
    double fund_freq = 2 * std::numbers::pi / sim.nbeads * sim.this_bead;
    if (sim.this_bead == 0) {
        pref = 1 / sqrt(sim.nbeads);
        std::fill(cart_to_nm_mat_row.begin(), cart_to_nm_mat_row.end(), pref);
    } else if (sim.this_bead < 0.5 * sim.nbeads) {
        pref = sqrt(2.0 / sim.nbeads);
        for (int i = 0; i < sim.nbeads; ++i)
            cart_to_nm_mat_row[i] = pref * cos(fund_freq * i);
    } else if (sim.this_bead == 0.5 * sim.nbeads) {
        pref = 1 / sqrt(sim.nbeads);
        for (int i = 0; i < sim.nbeads; ++i)
            cart_to_nm_mat_row[i] = pref * (i % 2 == 0 ? 1.0 : -1.0);
    } else {
        pref = sqrt(2.0 / sim.nbeads);
        for (int i = 0; i < sim.nbeads; ++i)
            cart_to_nm_mat_row[i] = -pref * sin(fund_freq * i);
    }
    
    // NM-to-Cartesian transformation matrix row
    pref = sqrt(2.0 / sim.nbeads);
    nm_to_cart_mat_row[0] = 1 / sqrt(sim.nbeads);
    for (int i = 1; i < 0.5 * sim.nbeads; ++i)
        nm_to_cart_mat_row[i] = pref * cos(fund_freq * i);
    if (sim.nbeads % 2 == 0)
        nm_to_cart_mat_row[sim.nbeads / 2] = 1 / sqrt(sim.nbeads) * (sim.this_bead % 2 == 0 ? 1.0 : -1.0);
    for (int i = std::ceil(0.5 * (sim.nbeads + 1)); i < sim.nbeads; ++i)
        nm_to_cart_mat_row[i] = -pref * sin(fund_freq * i);
}

NormalModes::~NormalModes() {
    MPI_Win_free(&win_coord_cartesian);
    MPI_Win_free(&win_coord_nm);
    MPI_Win_free(&win_momenta_cartesian);
    MPI_Win_free(&win_momenta_nm);
}

void NormalModes::shareData() {
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = globIndexAtom(axis, ptcl_idx);
            arr_coord_cartesian[glob_idx + sim.this_bead] = sim.coord(ptcl_idx, axis);
            arr_momenta_cartesian[glob_idx + sim.this_bead] = sim.momenta(ptcl_idx, axis);
        }
    }
}

double NormalModes::coordCarToNM(const int glob_idx) {
    double coord_nm = 0;
    for (int bead_idx = 0; bead_idx < sim.nbeads; ++bead_idx) {
        coord_nm += cart_to_nm_mat_row[bead_idx] * arr_coord_cartesian[glob_idx + bead_idx];
    }
    return coord_nm;
}

double NormalModes::momentumCarToNM(const int glob_idx) {
    double momentum_nm = 0;
    for (int bead_idx = 0; bead_idx < sim.nbeads; ++bead_idx) {
        momentum_nm += cart_to_nm_mat_row[bead_idx] * arr_momenta_cartesian[glob_idx + bead_idx];
    }
    return momentum_nm;
}

double NormalModes::coordNMToCar(const int glob_idx) {
    double coord_cartesian = 0;
    for (int bead_idx = 0; bead_idx < sim.nbeads; ++bead_idx) {
        coord_cartesian += nm_to_cart_mat_row[bead_idx] * arr_coord_nm[glob_idx + bead_idx];
    }
    return coord_cartesian;
}

double NormalModes::momentumNMToCar(const int glob_idx) {
    double momentum_cartesian = 0;
    for (int bead_idx = 0; bead_idx < sim.nbeads; ++bead_idx) {
        momentum_cartesian += nm_to_cart_mat_row[bead_idx] * arr_momenta_nm[glob_idx + bead_idx];
    }
    return momentum_cartesian;
}

void NormalModes::updateCartesianMomenta() {
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = globIndexAtom(axis, ptcl_idx);
            arr_momenta_cartesian[glob_idx + sim.this_bead] = momentumNMToCar(glob_idx);
            sim.momenta(ptcl_idx, axis) = arr_momenta_cartesian[glob_idx + sim.this_bead]; 
        }
    }
}

