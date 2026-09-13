#pragma once

#include "core/simulation_config.h"
#include "common.h"
#include "inireader.h"

#include <memory>

namespace Sections {
    constexpr const char* SYSTEM = "system";
    constexpr const char* SIMULATION = "simulation";
    constexpr const char* ACTION = "action";
    constexpr const char* EXT_POTENTIAL = "external_potential";
    constexpr const char* INT_POTENTIAL = "interaction_potential";
    constexpr const char* DUMP = "dump";
    constexpr const char* OBSERVABLES = "observables";
}

class Params {
public:
    explicit Params(const std::string& filename, int rank);
    [[nodiscard]] std::shared_ptr<SimulationConfig> load() const;
    double loadQuantity(
        const std::string& family,
        const std::string& section,
        const std::string& key,
        const std::string& default_value,
        double& destination,
        StringMap& units_destination
    ) const;
private:
    INIReader m_reader;
    int m_rank;

    void loadSimulationParams(SimulationConfig& config) const;
    void loadSystemParams(SimulationConfig& config) const;
    void loadPropagatorParams(SimulationConfig& config) const;
    void loadActionParams(SimulationConfig& config) const;
    void loadThermostatParams(SimulationConfig& config) const;
    void loadCoordInitParams(SimulationConfig& config) const;
    void loadVelocityInitParams(SimulationConfig& config) const;
    void loadExternalPotentialParams(SimulationConfig& config) const;
    void loadInteractionPotentialParams(SimulationConfig& config) const;
    void loadOutputParams(SimulationConfig& config) const;
    void loadObservableParams(SimulationConfig& config) const;
    void loadRpmdParams(SimulationConfig& config) const;
};
