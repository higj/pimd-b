#include "propagators/velocity_verlet.h"
#include "simulation.h"

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