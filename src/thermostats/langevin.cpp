#include "thermostats/langevin.h"
#include "thermostats/thermostat_coupling.h"
#include "simulation.h"
#include "common.h"
#include "normal_modes.h"

#include <numbers>

// A class for a Langevin thermostat 
LangevinThermostat::LangevinThermostat(Coupling& coupling, Params& param_obj, int rank) : Thermostat(coupling, param_obj) {
    double gamma;
    getVariant(param_obj.sim["gamma"], gamma);
    friction_coefficient = exp(-0.5 * gamma * dt);
    noise_coefficient = sqrt((1 - friction_coefficient * friction_coefficient) * mass / thermo_beta);
    int params_seed;
    getVariant(param_obj.sim["seed"], params_seed);
    mars_gen = std::make_unique<RanMars>(params_seed + rank);
}

void LangevinThermostat::momentaUpdate() {
    //std::normal_distribution<double> normal; // At the moment we use the Marsaglia generator

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double noise = mars_gen->gaussian(); // LAMMPS pimd/langevin random numbers
            double momentum_for_calc = coupling.getMomentumForCalc(ptcl_idx, axis);
            double& momentum_for_update = coupling.getMomentumForUpdate(ptcl_idx, axis);
            // Perturb the momenta with a Langevin thermostat
            momentum_for_update = friction_coefficient * momentum_for_calc + noise_coefficient * noise;
        }
    }
}