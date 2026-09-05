#pragma once

#include "distinguishable_primitive_kinetic_energy_strategy.h"
#include "contexts/bead_context.h"

#include <memory>

class BosonicExchangeBase;

class BosonicPrimitiveKineticEnergyStrategy : public DistinguishablePrimitiveKineticEnergyStrategy {
public:
    BosonicPrimitiveKineticEnergyStrategy(
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange
    );

    double calculateSpringContribution(
        const VecArray& coord,
        const VecArray& prev_coord,
        const SpringContext& spring_ctx, 
        const BoxContext& box_ctx,
        const BeadContext& bead_ctx
    ) override;

private:
    std::shared_ptr<BosonicExchangeBase> m_bosonic_exchange;
};
