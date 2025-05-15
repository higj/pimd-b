#include "propagators/velocity_verlet.h"
#include "simulation.h"

VelocityVerletPropagator::VelocityVerletPropagator(Params& param_obj, dVec& coord, dVec& momenta, dVec& forces) : 
    Propagator(param_obj, coord, momenta, forces) {
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
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            momenta(ptcl_idx, axis) += 0.5 * dt * forces(ptcl_idx, axis);
        }
    }
}

void VelocityVerletPropagator::coordsStep() {
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            coord(ptcl_idx, axis) += dt * momenta(ptcl_idx, axis) / mass;
        }
    }
}