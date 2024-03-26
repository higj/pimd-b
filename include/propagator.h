#pragma once

#include <vector>
#include "mpi.h"
#include "common.h"

class Simulation;

class Propagator {
public:
    explicit Propagator(Simulation& _sim);
    ~Propagator() = default;
    
    virtual void step() = 0;
protected:
    Simulation& sim; // Reference to the simulation object
};

/* -------------------------------- */

class VelocityVerletPropagator : public Propagator {
public:
    VelocityVerletPropagator(Simulation& _sim);
    
    void step() override;
};

/* -------------------------------- */

class NormalModesPropagator : public Propagator {
public:
    NormalModesPropagator(Simulation& _sim);
    ~NormalModesPropagator();
    void step() override;

private:
    int axis_stride, atom_stride;  // For indexing purposes
    double *arr_coord_cartesian, *arr_coord_nm, *arr_momenta_cartesian, *arr_momenta_nm;  // Arrays that contain the coordinates and momenta of the WHOLE system in two representations
    MPI_Win win_coord_cartesian, win_coord_nm, win_momenta_cartesian, win_momenta_nm;     // Window objects associated with the arrays
    double freq, c, s, m_omega;
    std::vector<double> cart_to_nm_mat_row, nm_to_cart_mat_row;
    dVec ext_forces, spring_forces;
    
    inline int _glob_idx(int axis, int atom, int bead);
    inline int _glob_idx_nobead(int axis, int atom);
};
