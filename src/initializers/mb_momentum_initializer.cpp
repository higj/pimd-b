#include "initializers/mb_momentum_initializer.h"
#include "core/random_generators.h"
//#include <random>

MaxwellBoltzmannMomentumInitializer::MaxwellBoltzmannMomentumInitializer(
    const std::shared_ptr<RandomGenerators>& rng,
    const std::shared_ptr<SystemState>& state,
    double mass,
    double thermo_beta) : MomentumInitializer(state, mass), m_rng(rng), m_thermo_beta(thermo_beta) {
}

/**
 * Initializes momenta according to the Maxwell-Boltzmann distribution.
 */
void MaxwellBoltzmannMomentumInitializer::initialize() {
    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            m_state->momenta(ptcl_idx, axis) = m_mass * sampleMaxwellBoltzmann(m_thermo_beta, m_mass);
        }
    }

    // When generating momenta from the Maxwell-Boltzmann distribution, zero the total momentum
    m_state->zeroMomentum();
}

/**
 * Samples velocities from the Maxwell-Boltzmann distribution.
 *
 * @return Sampled velocity from the Maxwell-Boltzmann distribution.
 */
double MaxwellBoltzmannMomentumInitializer::sampleMaxwellBoltzmann(double thermo_beta, double mass) const
{
    //std::normal_distribution<double> normal(0.0, 1 / sqrt(thermo_beta * mass));
    //return normal(rand_gen);
    return m_rng->normal(0.0, 1.0 / sqrt(thermo_beta * mass));
}