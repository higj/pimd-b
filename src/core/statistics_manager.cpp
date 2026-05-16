#include "core/statistics_manager.h"
#include "common.h"

#if FACTORIAL_BOSONIC_ALGORITHM
#include "bosonic_exchange/factorial_bosonic_exchange.h"
#else
#include "bosonic_exchange/quadratic_bosonic_exchange.h"
#endif

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

StatisticsManager& StatisticsManager::getInstance()
{
    static StatisticsManager instance;
    return instance;
}

void StatisticsManager::initializeBosonic(
    bool is_bosonic,
    const BeadContext& bead_ctx,
    const ThermalContext& thermal_ctx,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    const std::shared_ptr<SystemState>& state
)
{
    if (m_initialized)
    {
        throw std::runtime_error("StatisticsManager has already been initialized.");
    }

    m_is_bosonic = is_bosonic;
    m_bead_ctx = bead_ctx;

    const bool is_bosonic_bead = is_bosonic &&
        (m_bead_ctx.this_bead == 0 || m_bead_ctx.this_bead == m_bead_ctx.nbeads - 1);

    if (is_bosonic_bead)
    {
        m_bosonic_exchange = createBosonicExchangeObject(thermal_ctx, spring_ctx, box_ctx, state);
    }
    else
    {
        m_bosonic_exchange = nullptr;
    }

    m_initialized = true;
}

void StatisticsManager::checkInitialized() const {
    if (!m_initialized) {
        throw std::runtime_error("StatisticsManager not initialized!");
    }
}

std::shared_ptr<BosonicExchangeBase> StatisticsManager::createBosonicExchangeObject(
    const ThermalContext& thermal_ctx,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    const std::shared_ptr<SystemState>& state
)
{
    std::shared_ptr<VecArray> x_first_bead;
    std::shared_ptr<VecArray> x_last_bead;

    if (state->currentBead() == 0)
    {
        // At the first imaginary time slice, the last ("P") slice is the previous one
        x_first_bead = std::shared_ptr<VecArray>(state, &state->coord);
        x_last_bead = std::shared_ptr<VecArray>(state, &state->prev_coord);
    }
    else
    {
        // At the last imaginary time slice ("P"), the first slice is the next one
        x_first_bead = std::shared_ptr<VecArray>(state, &state->next_coord);
        x_last_bead = std::shared_ptr<VecArray>(state, &state->coord);
    }

#if FACTORIAL_BOSONIC_ALGORITHM
    return std::make_shared<FactorialBosonicExchange>(
        x_first_bead,
        x_last_bead,
        thermal_ctx,
        spring_ctx,
        box_ctx,
        m_bead_ctx
    );
#else
    return std::make_shared<BosonicExchange>(
        x_first_bead,
        x_last_bead,
        thermal_ctx,
        spring_ctx,
        box_ctx,
        m_bead_ctx
    );
#endif
}

std::unique_ptr<SpringForceStrategy> StatisticsManager::createSpringForceStrategy()
{
    checkInitialized();

    if (m_is_bosonic)
    {
        return std::make_unique<BosonicSpringForceStrategy>(m_bosonic_exchange, m_bead_ctx);
    }

    return std::make_unique<DistinguishableSpringForceStrategy>();
}

std::unique_ptr<PrimitiveKineticEnergyStrategy> StatisticsManager::createPrimitiveKineticEnergyStrategy()
{
    checkInitialized();

    if (m_bead_ctx.this_bead == 0 && m_is_bosonic)
    {
        return std::make_unique<BosonicPrimitiveKineticEnergyStrategy>(m_bosonic_exchange);
    }

    return std::make_unique<DistinguishablePrimitiveKineticEnergyStrategy>();
}

std::unique_ptr<ClassicalSpringEnergyStrategy> StatisticsManager::createClassicalSpringEnergyStrategy()
{
    checkInitialized();

    if (m_bead_ctx.this_bead == 0 && m_is_bosonic)
    {
        return std::make_unique<BosonicClassicalSpringEnergyStrategy>(m_bosonic_exchange);
    }

    return std::make_unique<DistinguishableClassicalSpringEnergyStrategy>();
}

std::unique_ptr<BosonicProbabilityStrategy> StatisticsManager::createBosonicProbabilityStrategy()
{
    checkInitialized();

    if (!m_is_bosonic)
    {
        throw std::runtime_error("BosonicProbabilityStrategy can only be used in bosonic simulations.");
    }

    if (m_bead_ctx.this_bead == 0)
    {
        return std::make_unique<ExteriorBosonicProbabilityStrategy>(m_bosonic_exchange);
    }

    return std::make_unique<InteriorBosonicProbabilityStrategy>();
}

std::unique_ptr<NormalModesMomentaStrategy> StatisticsManager::createNormalModesMomentaStrategy()
{
    checkInitialized();

    if (m_is_bosonic && (m_bead_ctx.this_bead == 0 || m_bead_ctx.this_bead == m_bead_ctx.nbeads - 1))
    {
        return std::make_unique<BosonicNormalModesMomentaStrategy>(m_bosonic_exchange);
    }

    return std::make_unique<DistinguishableNormalModesMomentaStrategy>();
}
