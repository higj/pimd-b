#include <numbers>
#include "common.h"
#include "normal_modes.h"
#include "params.h"

NormalModes::NormalModes(Params& param_obj, int this_bead, dVec& coord, dVec& momenta) :
    this_bead(this_bead), coord(coord), momenta(momenta)
{
    getVariant(param_obj.sys["natoms"], natoms);
    getVariant(param_obj.sim["nbeads"], nbeads);
    axis_stride = natoms*nbeads;
    atom_stride = nbeads;
    cart_to_nm_mat_row.resize(nbeads);
    nm_to_cart_mat_row.resize(nbeads);
    double temperature;
    getVariant(param_obj.sys["temperature"], temperature);
#if IPI_CONVENTION
    // i-Pi convention [J. Chem. Phys. 133, 124104 (2010)]
    double omega_p = nbeads * Constants::kB * temperature / Constants::hbar;
#else
    // Tuckerman convention
    double omega_p = sqrt(nbeads) * Constants::kB * temperature / Constants::hbar;
#endif

    // Allocate shared memory
    if (this_bead == 0) {
        MPI_Win_allocate_shared(natoms*nbeads*NDIM*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_coord_cartesian, &win_coord_cartesian);
        MPI_Win_allocate_shared(natoms*nbeads*NDIM*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_coord_nm, &win_coord_nm);
        MPI_Win_allocate_shared(natoms*nbeads*NDIM*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_momenta_cartesian, &win_momenta_cartesian);
        MPI_Win_allocate_shared(natoms*nbeads*NDIM*sizeof(double), sizeof(double),
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
    double freq = 2 * omega_p * sin(this_bead * std::numbers::pi / nbeads);
    
    // Cartesian-to-nm transformation matrix (one row because parallelized)
    double pref;
    double fund_freq = 2 * std::numbers::pi / nbeads * this_bead;
    if (this_bead == 0) {
        pref = 1 / sqrt(nbeads);
        std::fill(cart_to_nm_mat_row.begin(), cart_to_nm_mat_row.end(), pref);
    } else if (this_bead < 0.5 * nbeads) {
        pref = sqrt(2.0 / nbeads);
        for (int i = 0; i < nbeads; ++i)
            cart_to_nm_mat_row[i] = pref * cos(fund_freq * i);
    } else if (this_bead == 0.5 * nbeads) {
        pref = 1 / sqrt(nbeads);
        for (int i = 0; i < nbeads; ++i)
            cart_to_nm_mat_row[i] = pref * (i % 2 == 0 ? 1.0 : -1.0);
    } else {
        pref = sqrt(2.0 / nbeads);
        for (int i = 0; i < nbeads; ++i)
            cart_to_nm_mat_row[i] = -pref * sin(fund_freq * i);
    }
    
    // NM-to-Cartesian transformation matrix row
    pref = sqrt(2.0 / nbeads);
    nm_to_cart_mat_row[0] = 1 / sqrt(nbeads);
    for (int i = 1; i < 0.5 * nbeads; ++i)
        nm_to_cart_mat_row[i] = pref * cos(fund_freq * i);
    if (nbeads % 2 == 0)
        nm_to_cart_mat_row[nbeads / 2] = 1 / sqrt(nbeads) * (this_bead % 2 == 0 ? 1.0 : -1.0);
    for (int i = std::ceil(0.5 * (nbeads + 1)); i < nbeads; ++i)
        nm_to_cart_mat_row[i] = -pref * sin(fund_freq * i);
}

NormalModes::~NormalModes() {
    MPI_Win_free(&win_coord_cartesian);
    MPI_Win_free(&win_coord_nm);
    MPI_Win_free(&win_momenta_cartesian);
    MPI_Win_free(&win_momenta_nm);
}

void NormalModes::shareData() {
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = globIndexAtom(axis, ptcl_idx);
            arr_coord_cartesian[glob_idx + this_bead] = coord(ptcl_idx, axis);
            arr_momenta_cartesian[glob_idx + this_bead] = momenta(ptcl_idx, axis);
        }
    }
}

double NormalModes::coordCarToNM(const int glob_idx) {
    double coord_nm = 0;
    for (int bead_idx = 0; bead_idx < nbeads; ++bead_idx) {
        coord_nm += cart_to_nm_mat_row[bead_idx] * arr_coord_cartesian[glob_idx + bead_idx];
    }
    return coord_nm;
}

double NormalModes::momentumCarToNM(const int glob_idx) {
    double momentum_nm = 0;
    for (int bead_idx = 0; bead_idx < nbeads; ++bead_idx) {
        momentum_nm += cart_to_nm_mat_row[bead_idx] * arr_momenta_cartesian[glob_idx + bead_idx];
    }
    return momentum_nm;
}

double NormalModes::coordNMToCar(const int glob_idx) {
    double coord_cartesian = 0;
    for (int bead_idx = 0; bead_idx < nbeads; ++bead_idx) {
        coord_cartesian += nm_to_cart_mat_row[bead_idx] * arr_coord_nm[glob_idx + bead_idx];
    }
    return coord_cartesian;
}

double NormalModes::momentumNMToCar(const int glob_idx) {
    double momentum_cartesian = 0;
    for (int bead_idx = 0; bead_idx < nbeads; ++bead_idx) {
        momentum_cartesian += nm_to_cart_mat_row[bead_idx] * arr_momenta_nm[glob_idx + bead_idx];
    }
    return momentum_cartesian;
}

void NormalModes::updateCartesianMomenta() {
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = globIndexAtom(axis, ptcl_idx);
            arr_momenta_cartesian[glob_idx + this_bead] = momentumNMToCar(glob_idx);
            momenta(ptcl_idx, axis) = arr_momenta_cartesian[glob_idx + this_bead]; 
        }
    }
}

