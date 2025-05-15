#include "propagators/velocity_verlet.h"
#include "simulation.h"

VelocityVerletPropagator::VelocityVerletPropagator(Simulation& _sim) : Propagator(_sim) {
}

void VelocityVerletPropagator::preForceStep() {
    // First step: momenta are propagated by half a step ("B" step)
    momentStep();

    // Second step: positions are propagated using the new momenta ("A" step)
    coordsStep();
}
void VelocityVerletPropagator::postForceStep() {
    // Final step: momenta are propagated once more ("B" step)
    momentStep();
}

void VelocityVerletPropagator::momentStep() {
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            sim.momenta(ptcl_idx, axis) += 0.5 * sim.dt * sim.forces(ptcl_idx, axis);
        }
    }
}

void VelocityVerletPropagator::coordsStep() {
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            sim.coord(ptcl_idx, axis) += sim.dt * sim.momenta(ptcl_idx, axis) / sim.mass;
        }
    }
}