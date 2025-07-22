#pragma once

#include "contexts/bead_context.h"

#include <memory>

class SpringForceStrategy;
class PrimitiveKineticEnergyStrategy;
class ClassicalSpringEnergyStrategy;
class BosonicProbabilityStrategy;
class NormalModesMomentaStrategy;

class BosonicExchangeBase;
class SystemState;

struct ThermalContext;
struct SpringContext;
struct BoxContext;

class StatisticsManager {
public:
    // Singleton access
    static StatisticsManager& getInstance();

    // Delete copy constructor and assignment operator
    StatisticsManager(const StatisticsManager&) = delete;
    StatisticsManager& operator=(const StatisticsManager&) = delete;

    // Initialize the exchange for the current bead context
    // This encapsulates all the bosonic logic that was in Simulation::initializeExchange
    void initializeBosonic(
        bool is_bosonic,
        const BeadContext& bead_ctx,
        const ThermalContext& thermal_ctx,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        const std::shared_ptr<SystemState>& state
    );

    // Strategy factory methods - now instance methods that use internal state
    std::unique_ptr<SpringForceStrategy> createSpringForceStrategy();
    std::unique_ptr<PrimitiveKineticEnergyStrategy> createPrimitiveKineticEnergyStrategy();
    std::unique_ptr<ClassicalSpringEnergyStrategy> createClassicalSpringEnergyStrategy();
    std::unique_ptr<BosonicProbabilityStrategy> createBosonicProbabilityStrategy();
    std::unique_ptr<NormalModesMomentaStrategy> createNormalModesMomentaStrategy();

    // Query methods to check if bosonic statistics are active
    [[nodiscard]] bool isBosonic() const { return m_is_bosonic; }
    [[nodiscard]] bool isBosonicExchangeActive() const { return m_bosonic_exchange != nullptr; }

private:
    StatisticsManager() = default;
    ~StatisticsManager() = default;

    // Internal method to create the appropriate bosonic exchange object
    std::shared_ptr<BosonicExchangeBase> createBosonicExchangeObject(
        const ThermalContext& thermal_ctx,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        const std::shared_ptr<SystemState>& state
    );

    // Internal state
    std::shared_ptr<BosonicExchangeBase> m_bosonic_exchange;
    BeadContext m_bead_ctx;
    bool m_is_bosonic = false;
};