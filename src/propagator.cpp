#include "propagator.h"
#include "simulation.h"


Propagator::Propagator(Simulation& _sim) :
    sim(_sim) {
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

#if RECENTER
    // Recentering should be attempted only in the case of periodic boundary conditions.
    // Also, with the current implementation, it can only work for distinguishable particles.
    if (sim.pbc && !sim.bosonic) {
        // If the initial bead has moved outside of the fundamental cell,
        // then rigidly translate the whole polymer such that
        // the polymer will start in the fundamental cell.

        // Vector that stores by how much the polymers should be translated (same as Delta(w)*L)
        dVec shift(sim.natoms);

        if (sim.this_bead == 0) {
            // In the distinguishable case, the outer loop is the loop over all the ring polymers.
            for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
                for (int axis = 0; axis < NDIM; ++axis) {
                    shift(ptcl_idx, axis) = std::nearbyint(coord(ptcl_idx, axis) / sim.size);
                    // Shift the initial bead back to the fundamental cell
                    coord(ptcl_idx, axis) -= sim.size * shift(ptcl_idx, axis);
                }
            }
        }

        // Broadcast the shift vector from process 0 (first bead) to all other processes
        const int shift_vec_size = shift.size();
        MPI_Bcast(shift.data(), shift_vec_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // For other time-slices, we receive info from the first bead and make a decision whether to move the polymer
        if (sim.this_bead != 0) {
            for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
                for (int axis = 0; axis < NDIM; ++axis) {
                    // Shift the current bead according to the shift of the initial bead
                    sim.coord(ptcl_idx, axis) -= sim.size * shift(ptcl_idx, axis);
                }
            }
        }
    }
#endif

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

NormalModesPropagator::NormalModesPropagator(Simulation& _sim) : Propagator(_sim) {
}

void NormalModesPropagator::step() {
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

#if RECENTER
    // Recentering should be attempted only in the case of periodic boundary conditions.
    // Also, with the current implementation, it can only work for distinguishable particles.
    if (sim.pbc && !sim.bosonic) {
        // If the initial bead has moved outside of the fundamental cell,
        // then rigidly translate the whole polymer such that
        // the polymer will start in the fundamental cell.

        // Vector that stores by how much the polymers should be translated (same as Delta(w)*L)
        dVec shift(sim.natoms);

        if (sim.this_bead == 0) {
            // In the distinguishable case, the outer loop is the loop over all the ring polymers.
            for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
                for (int axis = 0; axis < NDIM; ++axis) {
                    shift(ptcl_idx, axis) = std::nearbyint(coord(ptcl_idx, axis) / sim.size);
                    // Shift the initial bead back to the fundamental cell
                    coord(ptcl_idx, axis) -= sim.size * shift(ptcl_idx, axis);
                }
            }
        }

        // Broadcast the shift vector from process 0 (first bead) to all other processes
        const int shift_vec_size = shift.size();
        MPI_Bcast(shift.data(), shift_vec_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // For other time-slices, we receive info from the first bead and make a decision whether to move the polymer
        if (sim.this_bead != 0) {
            for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
                for (int axis = 0; axis < NDIM; ++axis) {
                    // Shift the current bead according to the shift of the initial bead
                    sim.coord(ptcl_idx, axis) -= sim.size * shift(ptcl_idx, axis);
                }
            }
        }
    }
#endif

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