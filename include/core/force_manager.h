#pragma once

#include "system_state.h"
#include "potentials/potential.h"
#include "contexts/box_context.h"

#include <memory>

struct SpringContext;
class SpringForceStrategy;

class ForceManager {
public:
    /**
     * @param ext_potential Pre-built external potential (ownership transferred).
     * @param int_potential Pre-built interaction potential (ownership transferred).
     * @param cutoff        Interaction cutoff radius in internal units. Negative means no cutoff.
     * @param box_ctx       Box context used to clamp the cutoff to L/2 when PBC is active.
     */
    ForceManager(
        std::unique_ptr<Potential> ext_potential,
        std::unique_ptr<Potential> int_potential,
        double cutoff,
        const BoxContext& box_ctx
    );

    ~ForceManager();

    /**
     * Updates the physical forces acting on the particles. This includes both the forces
     * due external potentials and the interaction forces between the particles.
     *
     * @param state Object representing the current state of the system, including forces acting on particles.
     * @param box_ctx Context for the box, which may include periodic boundary conditions and other box-related parameters.
     */
    void updatePhysicalForces(SystemState& state, const BoxContext& box_ctx) const;

    /**
     * Updates both the spring and physical forces' arrays.
     *
     * @param state Object representing the current state of the system, including forces acting on particles.
     * @param spring_ctx Context containing parameters related to the springs, such as the spring constant.
     * @param box_ctx Context for the box, which may include periodic boundary conditions and other box-related parameters.
     */
    void updateForces(
        SystemState& state, 
        const SpringContext& spring_ctx, 
        const BoxContext& box_ctx
    ) const;

    std::unique_ptr<Potential> ext_potential;
    std::unique_ptr<Potential> int_potential;
    double cutoff = -1.0; // Interaction potential cutoff

private:
    std::unique_ptr<SpringForceStrategy> m_spring_force_strategy;
};