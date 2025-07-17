#pragma once

#include <memory>

#include "contexts/bead_context.h"

class BosonicExchangeBase;
class PrimitiveKineticEnergyStrategy;
class ClassicalSpringEnergyStrategy;
class SpringForceStrategy;
class BosonicProbabilityStrategy;
class NormalModesMomentaStrategy;

class StatisticsManager
{
public:
    static std::unique_ptr<SpringForceStrategy> createSpringForceStrategy(
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
        const BeadContext& bead_ctx,
        bool is_bosonic
    );

    static std::unique_ptr<PrimitiveKineticEnergyStrategy> createPrimitiveKineticEnergyStrategy(
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
        const BeadContext& bead_ctx,
        bool is_bosonic
    );

    static std::unique_ptr<ClassicalSpringEnergyStrategy> createClassicalSpringEnergyStrategy(
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
        const BeadContext& bead_ctx,
        bool is_bosonic
    );

    static std::unique_ptr<BosonicProbabilityStrategy> createBosonicProbabilityStrategy(
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
        const BeadContext& bead_ctx,
        bool is_bosonic
    );

    static std::unique_ptr<NormalModesMomentaStrategy> createNormalModesMomentaStrategy(
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
        const BeadContext& bead_ctx,
        bool is_bosonic
    );
};
