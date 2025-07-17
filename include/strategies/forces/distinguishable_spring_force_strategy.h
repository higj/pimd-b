#pragma once

#include "spring_force_strategy.h"

class DistinguishableSpringForceStrategy : public SpringForceStrategy {
public:
    /**
     * Updates the spring forces exerted on the beads of distinguishable particles,
     * which are given by Eqn. (12.6.4) in Tuckerman (1st ed).
     *
     * @param state Current state of the system, including coordinates and forces.
     * @param spring_ctx Context containing parameters related to the springs, such as the spring constant.
     * @param box_ctx Context containing parameters related to the simulation box, such as periodic boundary conditions.
     */
    void updateSpringForces(SystemState& state, const SpringContext& spring_ctx, const BoxContext& box_ctx) override;
};
