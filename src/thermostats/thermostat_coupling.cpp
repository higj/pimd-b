#include "thermostats/thermostat_coupling.h"
#include "propagators/normal_modes/normal_modes.h"

Coupling::Coupling(const std::shared_ptr<dVec>& momenta) : m_momenta(momenta) {}

/* -------------------------------- */

CartesianCoupling::CartesianCoupling(const std::shared_ptr<dVec>& momenta) : Coupling(momenta) {}

void CartesianCoupling::mpiCommunication() {}

double CartesianCoupling::getMomentumForCalc(const int ptcl_idx, const int axis) {
    return (*m_momenta)(ptcl_idx, axis);
}

double& CartesianCoupling::getMomentumForUpdate(const int ptcl_idx, const int axis) {
    return (*m_momenta)(ptcl_idx, axis);
} 

void CartesianCoupling::updateCoupledMomenta() {}

/* -------------------------------- */

NormalModesCoupling::NormalModesCoupling(const std::shared_ptr<dVec>& momenta, const std::shared_ptr<NormalModes>& normal_modes, int this_bead)
: Coupling(momenta), m_normal_modes(normal_modes), m_this_bead(this_bead) {}

void NormalModesCoupling::mpiCommunication() {
    m_normal_modes->shareData();
    MPI_Barrier(MPI_COMM_WORLD);
}
 
double NormalModesCoupling::getMomentumForCalc(const int ptcl_idx, const int axis) {
    const int glob_idx = m_normal_modes->globIndexAtom(axis, ptcl_idx);
    return m_normal_modes->coordCartesianToNormal(glob_idx);
} 
 
double& NormalModesCoupling::getMomentumForUpdate(const int ptcl_idx, const int axis) {
    const int glob_idx = m_normal_modes->globIndexAtom(axis, ptcl_idx);
    return m_normal_modes->arr_momenta_nm[glob_idx + m_this_bead];
}

void NormalModesCoupling::updateCoupledMomenta() {
    MPI_Barrier(MPI_COMM_WORLD);
    m_normal_modes->updateCartesianMomenta();
}
