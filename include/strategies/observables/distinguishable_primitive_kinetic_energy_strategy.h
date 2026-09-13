#pragma once

#include "primitive_kinetic_energy_strategy.h"

class DistinguishablePrimitiveKineticEnergyStrategy : public PrimitiveKineticEnergyStrategy {
public:
    double calculateSpringContribution(
        const VecArray& coord,
        const VecArray& prev_coord,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        const BeadContext& bead_ctx
    ) override;
};
