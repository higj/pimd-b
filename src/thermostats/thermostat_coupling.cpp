#include <memory>
#include "thermostats/thermostat_coupling.h"
#include "mpi.h"
#include "normal_modes.h"
#include "simulation.h"

Coupling::Coupling(dVec& _momenta) : momenta(_momenta) {}

/* -------------------------------- */

CartesianCoupling::CartesianCoupling(dVec& _momenta) : Coupling(_momenta) {}

void CartesianCoupling::mpiCommunication() {}

double CartesianCoupling::getMomentumForCalc(const int ptcl_idx, const int axis) {
    return momenta(ptcl_idx, axis);
}

double& CartesianCoupling::getMomentumForUpdate(const int ptcl_idx, const int axis) {
    return momenta(ptcl_idx, axis);
} 

void CartesianCoupling::updateCoupledMomenta() {}

/* -------------------------------- */

NMCoupling::NMCoupling(dVec& _momenta, NormalModes& normal_modes, int this_bead) : 
    Coupling(_momenta), normal_modes(normal_modes), this_bead(this_bead) {}

void NMCoupling::mpiCommunication() {
    normal_modes.shareData();
    MPI_Barrier(MPI_COMM_WORLD);
}
 
double NMCoupling::getMomentumForCalc(const int ptcl_idx, const int axis) {
    const int glob_idx = normal_modes.globIndexAtom(axis, ptcl_idx);
    return normal_modes.momentumCarToNM(glob_idx);
} 
 
double& NMCoupling::getMomentumForUpdate(const int ptcl_idx, const int axis) {
    const int glob_idx = normal_modes.globIndexAtom(axis, ptcl_idx);
    return normal_modes.arr_momenta_nm[glob_idx + this_bead];
}

void NMCoupling::updateCoupledMomenta() {
    MPI_Barrier(MPI_COMM_WORLD);
    normal_modes.updateCartesianMomenta();
}
