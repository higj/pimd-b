#include "core/system_state.h"
#include "mpi.h"

SystemState::SystemState(int rank, int nproc, int natoms, int nbeads) : m_rank(rank), m_nproc(nproc), m_natoms(natoms), m_nbeads(nbeads) {
    // Initialize the coordinate, momenta, and force arrays
    coord = dVec(natoms);
    prev_coord = dVec(natoms);
    next_coord = dVec(natoms);
    momenta = dVec(natoms);
    spring_forces = dVec(natoms);
    physical_forces = dVec(natoms);
}

void SystemState::zeroMomentum() {
    dVec momentum_cm;           // Resulting center of mass momentum vector
    dVec momentum_cm_per_bead;  // Contribution of the current time-slice to the center of mass momentum vector

    //const int natoms = momenta.len();

    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            momentum_cm_per_bead(0, axis) += momenta(ptcl_idx, axis);
        }
    }

    for (int axis = 0; axis < NDIM; ++axis) {
        momentum_cm_per_bead(0, axis) /= (m_natoms * m_nbeads);
    }

    MPI_Allreduce(
        momentum_cm_per_bead.data(), 
        momentum_cm.data(), 
        momentum_cm.size(), 
        MPI_DOUBLE, 
        MPI_SUM,
        MPI_COMM_WORLD
    );

    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            momenta(ptcl_idx, axis) -= momentum_cm(0, axis);
        }
    }
}

void SystemState::updatePreviousCoordinates(dVec& prev) {
    const int coord_size = coord.size();

    MPI_Sendrecv(
        coord.data(),
        coord_size,
        MPI_DOUBLE,
        (m_rank + 1) % m_nproc,
        0,
        prev.data(),
        coord_size,
        MPI_DOUBLE,
        (m_rank - 1 + m_nproc) % m_nproc,
        0,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE
    );

    // Ensure all processes have completed the neighbor communication
    MPI_Barrier(MPI_COMM_WORLD);
}

void SystemState::updateNextCoordinates(dVec& next) {
    const int coord_size = coord.size();

    MPI_Sendrecv(
        coord.data(),
        coord_size,
        MPI_DOUBLE,
        (m_rank - 1 + m_nproc) % m_nproc,
        0,
        next.data(),
        coord_size,
        MPI_DOUBLE,
        (m_rank + 1) % m_nproc,
        0,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE
    );

    // Ensure all processes have completed the neighbor communication
    MPI_Barrier(MPI_COMM_WORLD);
}

void SystemState::updateNeighboringCoordinates() {
    updatePreviousCoordinates(prev_coord);
    updateNextCoordinates(next_coord);
}

double SystemState::getTotalForce(int ptcl_idx, int axis) const {
#if IPI_CONVENTION
    // i-Pi convention; exp[-(beta/P)*H_cl]
    return spring_forces(ptcl_idx, axis) + physical_forces(ptcl_idx, axis);
#else
    // Corresponds to Eqn. (12.6.4) from Tuckerman; exp[-beta*H_cl]
    return spring_forces(ptcl_idx, axis) + physical_forces(ptcl_idx, axis) / m_nbeads;
#endif
}