#include "strategies/observables/bosonic_primitive_kinetic_energy_strategy.h"
#include "bosonic_exchange/bosonic_exchange_base.h"

BosonicPrimitiveKineticEnergyStrategy::BosonicPrimitiveKineticEnergyStrategy(
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange
) : m_bosonic_exchange(bosonic_exchange)
{
}

double BosonicPrimitiveKineticEnergyStrategy::calculateSpringContribution(
    const dVec& /* coord */,
    const dVec& /* prev_coord */,
    const SpringContext& /* spring_ctx */,
    const BoxContext& /* box_ctx */,
    const BeadContext& /* bead_ctx */
)
{
    return m_bosonic_exchange->primitiveEnergyEstimator();
}
