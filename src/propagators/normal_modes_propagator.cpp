#include "propagators/normal_modes_propagator.h"
#include "normal_modes.h"
#include "simulation.h"

#include <numbers>

NormalModesPropagator::NormalModesPropagator(Simulation& _sim, Params& param_obj, dVec& coord, dVec& momenta, dVec& forces,
                                            dVec& physical_forces, dVec& spring_forces) : 
    Propagator(_sim, param_obj, coord, momenta, forces),
    physical_forces(physical_forces),
    spring_forces(spring_forces)
{
    // Frequencies
    freq = 2 * sim.omega_p * sin(sim.this_bead * std::numbers::pi / sim.nbeads); 
    c = cos(freq * sim.dt);
    s = sin(freq * sim.dt);
    m_omega = sim.mass * freq;    
}

void NormalModesPropagator::preForceStep() {
    // Propagate momenta under external forces
    momentaExternalForces();

    sim.normal_modes->shareData();

    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            const int glob_idx = sim.normal_modes->globIndexAtom(axis, ptcl_idx);
            // Cartesian-to-nm transformation
            double coord_nm = sim.normal_modes->coordCarToNM(glob_idx);
            double momentum_nm = sim.normal_modes->momentumCarToNM(glob_idx);
            // Time propagation
            if (freq == 0) {
                sim.normal_modes->arr_coord_nm[glob_idx + sim.this_bead] = coord_nm + sim.dt / sim.mass * momentum_nm;
                sim.normal_modes->arr_momenta_nm[glob_idx + sim.this_bead] = momentum_nm;
            } else {
                sim.normal_modes->arr_coord_nm[glob_idx + sim.this_bead] = c * coord_nm + s / m_omega * momentum_nm;
                sim.normal_modes->arr_momenta_nm[glob_idx + sim.this_bead] = (-1) * m_omega * s * coord_nm + c * momentum_nm;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = sim.normal_modes->globIndexAtom(axis, ptcl_idx);
            // NM-to-Cartesian transformation
            double coord_cartesian = sim.normal_modes->coordNMToCar(glob_idx); 
            double momentum_cartesian = sim.normal_modes->momentumNMToCar(glob_idx);
            sim.normal_modes->arr_coord_cartesian[glob_idx + sim.this_bead] = coord_cartesian;
            sim.normal_modes->arr_momenta_cartesian[glob_idx + sim.this_bead] = momentum_cartesian;
        }
    }

    // Update forces
    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = sim.normal_modes->globIndexAtom(axis, ptcl_idx);
            sim.coord(ptcl_idx, axis) = sim.normal_modes->arr_coord_cartesian[glob_idx + sim.this_bead];
            sim.momenta(ptcl_idx, axis) = sim.normal_modes->arr_momenta_cartesian[glob_idx + sim.this_bead];
        }
    }
}
void NormalModesPropagator::postForceStep() {
    // Propagate momenta under external forces
    momentaExternalForces();
}

void NormalModesPropagator::momentaExternalForces() {
    if (!sim.bosonic) {
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
#if IPI_CONVENTION
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * physical_forces(ptcl_idx, axis);
#else
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * physical_forces(ptcl_idx, axis) / sim.nbeads;
#endif
    } else if (sim.this_bead == 0 || sim.this_bead == sim.nbeads - 1) {
        double inner_springs;
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                inner_springs = -sim.spring_constant * (2 * sim.coord(ptcl_idx, axis) - sim.prev_coord(ptcl_idx, axis) - sim.next_coord(ptcl_idx, axis));
#if IPI_CONVENTION
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * (physical_forces(ptcl_idx, axis) + spring_forces(ptcl_idx, axis) - inner_springs);
#else
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * (physical_forces(ptcl_idx, axis) / sim.nbeads + spring_forces(ptcl_idx, axis) - inner_springs);
#endif
            }
        }
    } else {
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
#if IPI_CONVENTION
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * physical_forces(ptcl_idx, axis);
#else
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * physical_forces(ptcl_idx, axis) / sim.nbeads;
#endif
    }
}
