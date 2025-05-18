#include "thermostats/langevin_thermostat.h"
#include "contexts/thermostats/thermostat_context.h"
#include "thermostats/thermostat_coupling.h"
#include "core/random_generators.h"

LangevinThermostat::LangevinThermostat(const ThermostatContext& context, const LangevinThermostatContext& langevin_context) :
Thermostat(context), m_langevin_context(langevin_context) {
    m_friction_coefficient = exp(-0.5 * langevin_context.gamma * context.dt);
    m_noise_coefficient = sqrt((1 - m_friction_coefficient * m_friction_coefficient) * context.mass / context.thermo_beta);
}

void LangevinThermostat::momentaUpdate() {
    //std::normal_distribution<double> normal; // At the moment we use the Marsaglia generator

    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double noise = m_langevin_context.rng->gaussian();
            double momentum_for_calc = coupling->getMomentumForCalc(ptcl_idx, axis);
            double& momentum_for_update = coupling->getMomentumForUpdate(ptcl_idx, axis);
            // Perturb the momenta with a Langevin thermostat
            momentum_for_update = m_friction_coefficient * momentum_for_calc + m_noise_coefficient * noise;
        }
    }
}