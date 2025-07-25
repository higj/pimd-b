#include "strategies/observables/bosonic_classical_spring_energy_strategy.h"
#include "bosonic_exchange/bosonic_exchange_base.h"

BosonicClassicalSpringEnergyStrategy::BosonicClassicalSpringEnergyStrategy(
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange
) : m_bosonic_exchange(bosonic_exchange)
{
}

double BosonicClassicalSpringEnergyStrategy::calculateSpringEnergy(
    const dVec& /* coord */,
    const dVec& /* prev_coord */,
    const SpringContext& /* spring_ctx */,
    const BoxContext& /* box_ctx */
)
{
    return m_bosonic_exchange->effectivePotential();
}
