#pragma once

#include "distinguishable_classical_spring_energy_strategy.h"

#include <memory>

class BosonicExchangeBase;

class BosonicClassicalSpringEnergyStrategy : public DistinguishableClassicalSpringEnergyStrategy {
public:
    BosonicClassicalSpringEnergyStrategy(
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange
    );

    double calculateSpringEnergy(
        const dVec& coord,
        const dVec& prev_coord,
        const SpringContext& spring_ctx, 
        const BoxContext& box_ctx
    ) override;

private:
    std::shared_ptr<BosonicExchangeBase> m_bosonic_exchange;
};
