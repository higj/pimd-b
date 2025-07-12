#include "thermostats/langevin_thermostat.h"
#include "thermostats/thermostat_coupling.h"
#include "core/random_generators.h"
#include "core/system_state.h"

LangevinThermostat::LangevinThermostat(
    const ThermalContext& thermal_ctx,
    const NormalModesContext& nm_ctx,
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<RandomGenerators>& rng,
    double gamma,
    double dt,
    double mass
) : Thermostat(thermal_ctx, nm_ctx, state), m_rng(rng), m_gamma(gamma), m_dt(dt), m_mass(mass), m_natoms(state->getNumAtoms()) {
    m_friction_coefficient = exp(-0.5 * gamma * m_dt);
    m_noise_coefficient = sqrt((1 - m_friction_coefficient * m_friction_coefficient) * m_mass / thermal_ctx.thermo_beta);
}

void LangevinThermostat::momentaUpdate() {
    //std::normal_distribution<double> normal; // At the moment we use the Marsaglia generator

    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double noise = m_rng->gaussian();
            double momentum_for_calc = m_coupling->getMomentumForCalc(ptcl_idx, axis);
            double& momentum_for_update = m_coupling->getMomentumForUpdate(ptcl_idx, axis);
            // Perturb the momenta with a Langevin thermostat
            momentum_for_update = m_friction_coefficient * momentum_for_calc + m_noise_coefficient * noise;
        }
    }
}