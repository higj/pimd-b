#include "strategies/observables/bosonic_primitive_kinetic_energy_strategy.h"
#include "bosonic_exchange/bosonic_exchange_base.h"

BosonicPrimitiveKineticEnergyStrategy::BosonicPrimitiveKineticEnergyStrategy(
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange
) : m_bosonic_exchange(bosonic_exchange)
{
}

double BosonicPrimitiveKineticEnergyStrategy::calculateSpringContribution(
    [[maybe_unused]] const dVec& coord,
    [[maybe_unused]] const dVec& prev_coord,
    [[maybe_unused]] const SpringContext& spring_ctx,
    [[maybe_unused]] const BoxContext& box_ctx,
    [[maybe_unused]] const BeadContext& bead_ctx
)
{
    return m_bosonic_exchange->primitiveEnergyEstimator();
}
