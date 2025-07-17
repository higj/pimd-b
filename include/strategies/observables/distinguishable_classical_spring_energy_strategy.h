#pragma once

#include "classical_spring_energy_strategy.h"

class DistinguishableClassicalSpringEnergyStrategy : public ClassicalSpringEnergyStrategy {
public:
    double calculateSpringEnergy(
        const dVec& coord,
        const dVec& prev_coord,
        const SpringContext& spring_ctx, 
        const BoxContext& box_ctx
    ) override;
};
