#pragma once

#include <memory>
#include "common.h"
#include "simulation_config.h"
#include "system_state.h"
#include "exchange_state.h"
#include "potentials/potential.h"

class ForceFieldManager {
public:
    ForceFieldManager(const std::shared_ptr<const SimulationConfig>& config);
    void updatePhysicalForces(SystemState& state) const;
    void updateSpringForces(SystemState& state, const ExchangeState& exchange_state) const;

    std::unique_ptr<Potential> m_ext_potential;
    std::unique_ptr<Potential> m_int_potential;
    double cutoff = -1.0; // Interaction potential cutoff

private:
    std::shared_ptr<const SimulationConfig> m_config;

    std::unique_ptr<Potential> initializePotential(const std::string& potential_name, const VariantMap& potential_options) const;

    /**
     * Helper function to apply minimum image if needed.
     */
    void applyMinimumImageIfNeeded(double& diff) const;

    /**
     * Helper function to calculate and accumulate spring force contribution
     */
    void addSpringForceContribution(SystemState& state, int ptcl_idx, int axis, double coord_diff) const;
    void updateDistinguishableSpringForces(SystemState& state) const;
    void updateBosonicSpringForces(SystemState& state, const ExchangeState& exchange_state) const;
};