#include <memory>
#include "thermostat_coupling.h"
#include "mpi.h"
#include "normal_modes.h"
#include "simulation.h"

Coupling::Coupling(Simulation& _sim) : sim(_sim) {}

/* -------------------------------- */

CartesianCoupling::CartesianCoupling(Simulation& _sim) : Coupling(_sim) {}

void CartesianCoupling::mpiCommunication() {}

double CartesianCoupling::getMomentumForCalc(const int ptcl_idx, const int axis) {
    return sim.momenta(ptcl_idx, axis);
}

double& CartesianCoupling::getMomentumForUpdate(const int ptcl_idx, const int axis) {
    return sim.momenta(ptcl_idx, axis);
} 

void CartesianCoupling::updateCoupledMomenta() {}

/* -------------------------------- */

NMCoupling::NMCoupling(Simulation& _sim) : Coupling(_sim) {}

void NMCoupling::mpiCommunication() {
    sim.normal_modes->shareData();
    MPI_Barrier(MPI_COMM_WORLD);
}
 
double NMCoupling::getMomentumForCalc(const int ptcl_idx, const int axis) {
    const int glob_idx = sim.normal_modes->globIndexAtom(axis, ptcl_idx);
    return sim.normal_modes->momentumCarToNM(glob_idx);
} 
 
double& NMCoupling::getMomentumForUpdate(const int ptcl_idx, const int axis) {
    const int glob_idx = sim.normal_modes->globIndexAtom(axis, ptcl_idx);
    return sim.normal_modes->arr_momenta_nm[glob_idx + sim.this_bead];
}

void NMCoupling::updateCoupledMomenta() {
    MPI_Barrier(MPI_COMM_WORLD);
    sim.normal_modes->updateCartesianMomenta();
}
