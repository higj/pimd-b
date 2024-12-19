#include "propagator.h"
#include <numbers>
#include "normal_modes.h"
#include "simulation.h"
#include "common.h"

Propagator::Propagator(Simulation& _sim) : sim(_sim) {
}

VelocityVerletPropagator::VelocityVerletPropagator(Simulation& _sim) : Propagator(_sim) {
}

void VelocityVerletPropagator::step() {
    // First step: momenta are propagated by half a step ("B" step)
    momentStep();

    // Second step: positions are propagated using the new momenta ("A" step)
    coordsStep();

    // Remember to update the neighboring coordinates after every coordinate propagation
    sim.updateNeighboringCoordinates();

    // Third step: forces are updated using the new positions
    sim.updateForces();

    // Fourth step: momenta are propagated once more ("B" step)
    momentStep();
}

void Propagator::momentStep() {
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * sim.forces(ptcl_idx, axis);
        }
    }
}

void Propagator::coordsStep() {
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            sim.coord(ptcl_idx, axis) += sim.dt * sim.momenta(ptcl_idx, axis) / sim.mass;
        }
    }
}

NormalModesPropagator::NormalModesPropagator(Simulation& _sim) : 
    Propagator(_sim),
    ext_forces(_sim.natoms),
    spring_forces(_sim.natoms)
{
    // Frequencies
    freq = 2 * sim.omega_p * sin(sim.this_bead * std::numbers::pi / sim.nbeads); 
    c = cos(freq * sim.dt);
    s = sin(freq * sim.dt);
    m_omega = sim.mass * freq;    
}

void NormalModesPropagator::step() {
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
    sim.updateNeighboringCoordinates();
    sim.updateForces();
    sim.updatePhysicalForces(ext_forces);
    sim.updateSpringForces(spring_forces);
    
    // Propagate momenta under external forces
    momentaExternalForces();
}

void NormalModesPropagator::momentaExternalForces() {
    if (!sim.bosonic) {
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
#if IPI_CONVENTION
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * ext_forces(ptcl_idx, axis);
#else
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * ext_forces(ptcl_idx, axis) / sim.nbeads;
#endif
    } else if (sim.this_bead == 0 || sim.this_bead == sim.nbeads - 1) {
        double inner_springs;
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                inner_springs = -sim.spring_constant * (2 * sim.coord(ptcl_idx, axis) - sim.prev_coord(ptcl_idx, axis) - sim.next_coord(ptcl_idx, axis));
#if IPI_CONVENTION
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * (ext_forces(ptcl_idx, axis) + spring_forces(ptcl_idx, axis) - inner_springs);
#else
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * (ext_forces(ptcl_idx, axis) / sim.nbeads + spring_forces(ptcl_idx, axis) - inner_springs);
#endif
            }
        }
    } else {
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
#if IPI_CONVENTION
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * ext_forces(ptcl_idx, axis);
#else
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * ext_forces(ptcl_idx, axis) / sim.nbeads;
#endif
    }
}
