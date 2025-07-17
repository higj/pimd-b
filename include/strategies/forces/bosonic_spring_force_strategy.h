#pragma once

#include "distinguishable_spring_force_strategy.h"
#include "contexts/bead_context.h"

#include <memory>

class BosonicExchangeBase;

class BosonicSpringForceStrategy : public DistinguishableSpringForceStrategy {
public:
    BosonicSpringForceStrategy(
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange, 
        const BeadContext& bead_ctx
    );

    void updateSpringForces(SystemState& state, const SpringContext& spring_ctx, const BoxContext& box_ctx) override;

private:
    /**
     * Updates the spring forces exerted on the beads of bosonic particles.
     * In the bosonic case, by default, the forces are evaluated using the algorithm
     * described in https://doi.org/10.1063/5.0173749. It is also possible to perform the
     * bosonic simulation using the original (inefficient) algorithm, that takes into
     * account all the N! permutations, by setting FACTORIAL_BOSONIC_ALGORITHM to true.
     *
     * @param state Object representing the current state of the system, including forces acting on particles.
     * @param spring_ctx Context containing parameters related to the springs, such as the spring constant.
     * @param box_ctx Context containing parameters related to the simulation box, such as periodic boundary conditions.
     */
    void updateBosonicForces(SystemState& state, const SpringContext& spring_ctx, const BoxContext& box_ctx);

    std::shared_ptr<BosonicExchangeBase> m_bosonic_exchange;
    BeadContext m_bead_ctx;

    bool m_is_exterior_bead;     // Is the current bead either 1 or P?
};
