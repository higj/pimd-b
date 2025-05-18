#pragma once

#include "core/simulation_config.h"
#include "common.h"
#include "inireader.h"

#include <memory>

namespace Sections {
    const std::string SYSTEM = "system";
    const std::string SIMULATION = "simulation";
    const std::string EXT_POTENTIAL = "external_potential";
    const std::string INT_POTENTIAL = "interaction_potential";
    const std::string DUMP = "dump";
    const std::string OBSERVABLES = "observables";
}

class Params {
public:
    explicit Params(const std::string& filename, const int& rank);
    [[nodiscard]] std::shared_ptr<SimulationConfig> load() const;
private:
    INIReader m_reader;
    int m_rank;

    void loadSimulationParams(SimulationConfig& config) const;
    void loadSystemParams(SimulationConfig& config) const;
    void loadPropagatorParams(SimulationConfig& config) const;
    void loadThermostatParams(SimulationConfig& config) const;
    void loadCoordInitParams(SimulationConfig& config) const;
    void loadMomentaInitParams(SimulationConfig& config) const;
    void loadExternalPotentialParams(SimulationConfig& config) const;
    void loadInteractionPotentialParams(SimulationConfig& config) const;
    void loadOutputParams(SimulationConfig& config) const;
    void loadObservableParams(SimulationConfig& config) const;
};
