#pragma once

#include <vector>
#include "mpi.h"
#include "common.h"

class Simulation;

class NormalModes {
public:
    NormalModes(Simulation& _sim);
    ~NormalModes();
    double *arr_coord_cartesian, *arr_coord_nm, *arr_momenta_cartesian, *arr_momenta_nm;  // Arrays that contain the coordinates and momenta of the WHOLE system in two representations
    std::vector<double> cart_to_nm_mat_row, nm_to_cart_mat_row;
    void shareData();
    double coordCarToNM(const int glob_idx);
    double momentumCarToNM(const int glob_idx);
    double coordNMToCar(const int glob_idx);
    double momentumNMToCar(const int glob_idx);
    int globIndexAtom(const int axis, int atom) const { return axis * axis_stride + atom * atom_stride; }
    void updateCartesianMomenta(); 
private:
    int axis_stride, atom_stride;  // For indexing purposes
    MPI_Win win_coord_cartesian, win_coord_nm, win_momenta_cartesian, win_momenta_nm;     // Window objects associated with the arrays
    int globIndexBead(int axis, int atom, int bead) const { return globIndexAtom(axis, atom) + bead; }
    Simulation& sim; // Reference to the simulation object
};
