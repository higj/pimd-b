#include "thermostats/langevin.h"
#include "thermostats/thermostat_coupling.h"
#include "simulation.h"
#include "common.h"
#include "normal_modes.h"

#include <numbers>

// A class for a Langevin thermostat 
LangevinThermostat::LangevinThermostat(Simulation& _sim, bool normal_modes) : Thermostat(_sim, normal_modes) {
    friction_coefficient = exp(-0.5 * sim.gamma * sim.dt);
    noise_coefficient = sqrt((1 - friction_coefficient * friction_coefficient) * sim.mass / sim.thermo_beta);
}

void LangevinThermostat::momentaUpdate() {
    //std::normal_distribution<double> normal; // At the moment we use the Marsaglia generator

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double noise = sim.mars_gen->gaussian(); // LAMMPS pimd/langevin random numbers
            double momentum_for_calc = coupling->getMomentumForCalc(ptcl_idx, axis);
            double& momentum_for_update = coupling->getMomentumForUpdate(ptcl_idx, axis);
            // Perturb the momenta with a Langevin thermostat
            momentum_for_update = friction_coefficient * momentum_for_calc + noise_coefficient * noise;
        }
    }
}