#include "propagator.h"

#include <numbers>

#include "simulation.h"
#include "common.h"

Propagator::Propagator(Simulation& _sim) : sim(_sim) {
}

VelocityVerletPropagator::VelocityVerletPropagator(Simulation& _sim) : Propagator(_sim) {
}

void VelocityVerletPropagator::step() {
    // First step: momenta are propagated by half a step ("B" step)
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * sim.forces(ptcl_idx, axis);
        }
    }

    // Second step: positions are propagated using the new momenta ("A" step)
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            sim.coord(ptcl_idx, axis) += sim.dt * sim.momenta(ptcl_idx, axis) / sim.mass;
        }
    }

    // Remember to update the neighboring coordinates after every coordinate propagation
    sim.updateNeighboringCoordinates();

    // Third step: forces are updated using the new positions
    sim.updateForces();

    // Fourth step: momenta are propagated once more ("B" step)
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * sim.forces(ptcl_idx, axis);
        }
    }
}

NormalModesPropagator::NormalModesPropagator(Simulation& _sim) : 
    Propagator(_sim),
    axis_stride(_sim.natoms*_sim.nbeads),
    atom_stride(_sim.nbeads),
    cart_to_nm_mat_row(_sim.nbeads),
    nm_to_cart_mat_row(_sim.nbeads),
    ext_forces(_sim.natoms),
    spring_forces(_sim.natoms)
{
    // Allocate shared memory
    if (sim.this_bead == 0) {
        MPI_Win_allocate_shared(sim.natoms*sim.nbeads*NDIM*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_coord_cartesian, &win_coord_cartesian);
        MPI_Win_allocate_shared(sim.natoms*sim.nbeads*NDIM*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_coord_nm, &win_coord_nm);
        MPI_Win_allocate_shared(sim.natoms*sim.nbeads*NDIM*sizeof(double), sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &arr_momenta_cartesian, &win_momenta_cartesian);
        MPI_Win_allocate_shared(sim.natoms*sim.nbeads*NDIM*sizeof(double), sizeof(double),
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
    freq = 2 * sim.omega_p * sin(sim.this_bead * std::numbers::pi / sim.nbeads);
    
    // Cartesian-to-nm transformation matrix (one row because parallelized)
    double pref;
    double fund_freq = 2 * std::numbers::pi / sim.nbeads * sim.this_bead;
    if (sim.this_bead == 0) {
        pref = 1 / sqrt(sim.nbeads);
        std::fill(cart_to_nm_mat_row.begin(), cart_to_nm_mat_row.end(), pref);
    } else if (sim.this_bead < 0.5 * sim.nbeads) {
        pref = sqrt(2.0 / sim.nbeads);
        for (int i = 0; i < sim.nbeads; ++i)
            cart_to_nm_mat_row[i] = pref * cos(fund_freq * i);
    } else if (sim.this_bead == 0.5 * sim.nbeads) {
        pref = 1 / sqrt(sim.nbeads);
        for (int i = 0; i < sim.nbeads; ++i)
            cart_to_nm_mat_row[i] = pref * (i % 2 == 0 ? 1.0 : -1.0);
    } else {
        pref = sqrt(2.0 / sim.nbeads);
        for (int i = 0; i < sim.nbeads; ++i)
            cart_to_nm_mat_row[i] = -pref * sin(fund_freq * i);
    }
    
    // NM-to-Cartesian transformation matrix row
    pref = sqrt(2.0 / sim.nbeads);
    nm_to_cart_mat_row[0] = 1 / sqrt(sim.nbeads);
    for (int i = 1; i < 0.5 * sim.nbeads; ++i)
        nm_to_cart_mat_row[i] = pref * cos(fund_freq * i);
    if (sim.nbeads % 2 == 0)
        nm_to_cart_mat_row[sim.nbeads / 2] = 1 / sqrt(sim.nbeads) * (sim.this_bead % 2 == 0 ? 1.0 : -1.0);
    for (int i = std::ceil(0.5 * (sim.nbeads + 1)); i < sim.nbeads; ++i)
        nm_to_cart_mat_row[i] = -pref * sin(fund_freq * i);
    
    c = cos(freq * sim.dt);
    s = sin(freq * sim.dt);
    m_omega = sim.mass * freq;    
}

NormalModesPropagator::~NormalModesPropagator() {
    MPI_Win_free(&win_coord_cartesian);
    MPI_Win_free(&win_coord_nm);
    MPI_Win_free(&win_momenta_cartesian);
    MPI_Win_free(&win_momenta_nm);
}

void NormalModesPropagator::step() {
    // Propagate momenta under external forces
    if (!sim.bosonic) {
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * ext_forces(ptcl_idx, axis);
    } else if (sim.this_bead == 0 || sim.this_bead == sim.nbeads - 1) {
        double inner_springs;
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                inner_springs = -sim.spring_constant * (2 * sim.coord(ptcl_idx, axis) - sim.prev_coord(ptcl_idx, axis) - sim.next_coord(ptcl_idx, axis));
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * (ext_forces(ptcl_idx, axis) + spring_forces(ptcl_idx, axis) - inner_springs);
            }
        }
    } else {
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * ext_forces(ptcl_idx, axis);
    }

    // Share data with other processes
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = globIndexAtom(axis, ptcl_idx);
            arr_coord_cartesian[glob_idx + sim.this_bead] = sim.coord(ptcl_idx, axis);
            arr_momenta_cartesian[glob_idx + sim.this_bead] = sim.momenta(ptcl_idx, axis);
        }
    }

    // Cartesian-to-nm transformation
    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = globIndexAtom(axis, ptcl_idx);
            // Cartesian-to-nm transformation
            double coord_nm = 0, momentum_nm = 0;
            for (int bead_idx = 0; bead_idx < sim.nbeads; ++bead_idx) {
                coord_nm += cart_to_nm_mat_row[bead_idx] * arr_coord_cartesian[glob_idx + bead_idx];
                momentum_nm += cart_to_nm_mat_row[bead_idx] * arr_momenta_cartesian[glob_idx + bead_idx];
            }
            // Time propagation
            if (freq == 0) {
                arr_coord_nm[glob_idx + sim.this_bead] = coord_nm + sim.dt / sim.mass * momentum_nm;
                arr_momenta_nm[glob_idx + sim.this_bead] = momentum_nm;
            } else {
                arr_coord_nm[glob_idx + sim.this_bead] = c * coord_nm + s / m_omega * momentum_nm;
                arr_momenta_nm[glob_idx + sim.this_bead] = (-1) * m_omega * s * coord_nm + c * momentum_nm;
            }
        }
    }

    // NM-to-Cartesian transformation
    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = globIndexAtom(axis, ptcl_idx);
            double coord_cartesian = 0, momentum_cartesian = 0;
            for (int bead_idx = 0; bead_idx < sim.nbeads; ++bead_idx) {
                coord_cartesian += nm_to_cart_mat_row[bead_idx] * arr_coord_nm[glob_idx + bead_idx];
                momentum_cartesian += nm_to_cart_mat_row[bead_idx] * arr_momenta_nm[glob_idx + bead_idx];
            }
            arr_coord_cartesian[glob_idx + sim.this_bead] = coord_cartesian;
            arr_momenta_cartesian[glob_idx + sim.this_bead] = momentum_cartesian;
        }
    }

    // Update forces
    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = globIndexAtom(axis, ptcl_idx);
            sim.coord(ptcl_idx, axis) = arr_coord_cartesian[glob_idx + sim.this_bead];
            sim.momenta(ptcl_idx, axis) = arr_momenta_cartesian[glob_idx + sim.this_bead];
        }
    }
    sim.updateNeighboringCoordinates();
    sim.updateForces();
    sim.updatePhysicalForces(ext_forces);
    sim.updateSpringForces(spring_forces);

    // Propagate momenta under external forces
    if (!sim.bosonic) {
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * ext_forces(ptcl_idx, axis);
    } else if (sim.this_bead == 0 || sim.this_bead == sim.nbeads - 1) {
        double inner_springs;
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                inner_springs = -sim.spring_constant * (2 * sim.coord(ptcl_idx, axis) - sim.prev_coord(ptcl_idx, axis) - sim.next_coord(ptcl_idx, axis));
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * (ext_forces(ptcl_idx, axis) + spring_forces(ptcl_idx, axis) - inner_springs);
            }
        }
    } else {
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
                sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * ext_forces(ptcl_idx, axis);
    }
}