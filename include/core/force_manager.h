#pragma once

#include <memory>
#include "common.h"
#include "simulation_config.h"
#include "system_state.h"
#include "potentials/potential.h"

class BosonicExchangeBase;
struct BeadContext;
struct BoxContext;
struct SpringContext;
class SpringForceStrategy;

class ForceManager {
public:
    explicit ForceManager(
        const SimulationConfig& config,
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
        const BeadContext& bead_ctx
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

    /**
     * Initializes the potential based on the input parameters.
     *
     * @param config Simulation configuration containing potential parameters.
     * @param potential_name Name of the potential.
     * @param potential_options Physical parameters of the potential.
     * @return Pointer to the initialized potential.
     */
    [[nodiscard]] static std::unique_ptr<Potential> initializePotential(
        const SimulationConfig& config,
        const std::string& potential_name,
        const VariantMap& potential_options
    );
};