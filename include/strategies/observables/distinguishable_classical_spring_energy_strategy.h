#pragma once

#include "classical_spring_energy_strategy.h"

class DistinguishableClassicalSpringEnergyStrategy : public ClassicalSpringEnergyStrategy {
public:
    double calculateSpringEnergy(
        const VecArray& coord,
        const VecArray& prev_coord,
        const SpringContext& spring_ctx, 
        const BoxContext& box_ctx
    ) override;
};
