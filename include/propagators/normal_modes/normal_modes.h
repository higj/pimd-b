#pragma once

#include <vector>
#include "mpi.h"

#include "contexts/normal_modes_context.h"

class NormalModes {
public:
    NormalModes(const NormalModesContext& context);
    ~NormalModes();
    double *arr_coord_cartesian, *arr_coord_nm, *arr_momenta_cartesian, *arr_momenta_nm;  // Arrays that contain the coordinates and momenta of the WHOLE system in two representations
    std::vector<double> cart_to_nm_mat_row, nm_to_cart_mat_row;
    
    void shareData();
    double coordCartesianToNormal(int glob_idx) const;
    double momentumCartesianToNormal(int glob_idx) const;
    double coordNormalToCartesian(int glob_idx) const;
    double momentumNormalToCartesian(int glob_idx) const;
    [[nodiscard]] int globIndexAtom(const int axis, int atom) const { return axis * axis_stride + atom * atom_stride; }
    void updateCartesianMomenta(); 
private:
    // Structure to group related MPI window and array pointers
    struct SharedMemory {
        double* array;
        MPI_Win window;
    };

    static void allocateSharedMemory(SharedMemory& mem, size_t size, int this_bead);
    void allocateAllSharedMemory(int this_bead);

    int axis_stride, atom_stride;  // For indexing purposes
    // Window objects associated with the arrays
    MPI_Win win_coord_cartesian, win_coord_nm, win_momenta_cartesian, win_momenta_nm;
    [[nodiscard]] int globIndexBead(int axis, int atom, int bead) const { return globIndexAtom(axis, atom) + bead; }
    
    NormalModesContext m_context;
};
