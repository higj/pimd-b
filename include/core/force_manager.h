#pragma once

#include <memory>
#include "common.h"
#include "simulation_config.h"
#include "system_state.h"
#include "exchange_state.h"
#include "potentials/potential.h"

class ForceManager {
public:
    explicit ForceManager(const std::shared_ptr<const SimulationConfig>& config);

    /**
     * Updates the physical forces acting on the particles. This includes both the forces
     * due external potentials and the interaction forces between the particles.
     *
     * @param state Object representing the current state of the system, including forces acting on particles.
     */
    void updatePhysicalForces(SystemState& state) const;

    /**
     * Updates the spring forces array.
     *
     * @param state Object representing the current state of the system, including forces acting on particles.
     * @param exchange_state Object representing the state of the exchange algorithm.
     */
    void updateSpringForces(
        SystemState& state, 
        const ExchangeState& exchange_state
    ) const;

    /**
     * Updates both the spring and physical forces' arrays.
     *
     * @param state Object representing the current state of the system, including forces acting on particles.
     * @param exchange_state Object representing the state of the exchange algorithm.
     */
    void updateForces(
        SystemState& state,
        const ExchangeState& exchange_state
    ) const;

    std::unique_ptr<Potential> ext_potential;
    std::unique_ptr<Potential> int_potential;
    double cutoff = -1.0; // Interaction potential cutoff

private:
    std::shared_ptr<const SimulationConfig> m_config;

    /**
     * Initializes the potential based on the input parameters.
     *
     * @param potential_name Name of the potential.
     * @param potential_options Physical parameters of the potential.
     * @return Pointer to the initialized potential.
     */
    [[nodiscard]] std::unique_ptr<Potential> initializePotential(
        const std::string& potential_name, 
        const VariantMap& potential_options
    ) const;

    /**
     * Helper function to apply minimum image if needed.
     */
    void applyMinimumImageIfNeeded(double& diff) const;

    /**
     * Helper function to calculate and accumulate spring force contribution
     */
    void addSpringForceContribution(
        SystemState& state, 
        int ptcl_idx, 
        int axis, 
        double coord_diff
    ) const;

    /**
     * Updates the spring forces exerted on the beads.
     * In the distinguishable case, the force is given by Eqn. (12.6.4) in Tuckerman (1st ed).
     *
     * @param state Object representing the current state of the system, including forces acting on particles.
     */
    void updateDistinguishableSpringForces(SystemState& state) const;

    /**
     * Updates the spring forces exerted on the beads.
     * In the bosonic case, by default, the forces are evaluated using the algorithm
     * described in https://doi.org/10.1063/5.0173749. It is also possible to perform the
     * bosonic simulation using the original (inefficient) algorithm, that takes into
     * account all the N! permutations, by setting FACTORIAL_BOSONIC_ALGORITHM to true.
     *
     * @param state Object representing the current state of the system, including forces acting on particles.
     * @param exchange_state Object representing the state of the exchange algorithm.
     */
    void updateBosonicSpringForces(
        SystemState& state,
        const ExchangeState& exchange_state
    ) const;
};