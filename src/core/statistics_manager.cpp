#include "core/statistics_manager.h"

#include "strategies/forces/distinguishable_spring_force_strategy.h"
#include "strategies/forces/bosonic_spring_force_strategy.h"

#include "strategies/observables/distinguishable_classical_spring_energy_strategy.h"
#include "strategies/observables/bosonic_classical_spring_energy_strategy.h"

#include "strategies/observables/distinguishable_primitive_kinetic_energy_strategy.h"
#include "strategies/observables/bosonic_primitive_kinetic_energy_strategy.h"

#include "strategies/observables/interior_bosonic_probability_strategy.h"
#include "strategies/observables/exterior_bosonic_probability_strategy.h"

#include "strategies/propagators/bosonic_normal_modes_momenta_strategy.h"
#include "strategies/propagators/distinguishable_normal_modes_momenta_strategy.h"

std::unique_ptr<SpringForceStrategy> StatisticsManager::createSpringForceStrategy(
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
    const BeadContext& bead_ctx,
    bool is_bosonic)
{
    if (is_bosonic)
    {
        return std::make_unique<BosonicSpringForceStrategy>(bosonic_exchange, bead_ctx);
    }

    return std::make_unique<DistinguishableSpringForceStrategy>();
}

std::unique_ptr<PrimitiveKineticEnergyStrategy> StatisticsManager::createPrimitiveKineticEnergyStrategy(
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange, 
    const BeadContext& bead_ctx,
    bool is_bosonic
) {
    if (bead_ctx.this_bead == 0 && is_bosonic) {
        return std::make_unique<BosonicPrimitiveKineticEnergyStrategy>(bosonic_exchange);
    }

    return std::make_unique<DistinguishablePrimitiveKineticEnergyStrategy>();
}

std::unique_ptr<ClassicalSpringEnergyStrategy> StatisticsManager::createClassicalSpringEnergyStrategy(
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
    const BeadContext& bead_ctx,
    bool is_bosonic)
{
    if (bead_ctx.this_bead == 0 && is_bosonic)
    {
        return std::make_unique<BosonicClassicalSpringEnergyStrategy>(bosonic_exchange);
    }

    return std::make_unique<DistinguishableClassicalSpringEnergyStrategy>();
}

std::unique_ptr<BosonicProbabilityStrategy> StatisticsManager::createBosonicProbabilityStrategy(
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
    const BeadContext& bead_ctx,
    bool is_bosonic)
{
    if (!is_bosonic) {
        throw std::runtime_error("BosonicProbabilityStrategy can only be used in bosonic simulations.");
    }

    if (bead_ctx.this_bead == 0)
    {
        return std::make_unique<ExteriorBosonicProbabilityStrategy>(bosonic_exchange);
    }

    return std::make_unique<InteriorBosonicProbabilityStrategy>();
}

std::unique_ptr<NormalModesMomentaStrategy> StatisticsManager::createNormalModesMomentaStrategy(
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
    const BeadContext& bead_ctx,
    bool is_bosonic)
{
    if (is_bosonic && (bead_ctx.this_bead == 0 || bead_ctx.this_bead == bead_ctx.nbeads - 1))
    {
        return std::make_unique<BosonicNormalModesMomentaStrategy>(bosonic_exchange);
    }

    return std::make_unique<DistinguishableNormalModesMomentaStrategy>();
}
