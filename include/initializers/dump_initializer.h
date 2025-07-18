#pragma once

#include "dumps.h"
#include "contexts/box_context.h"
#include "contexts/velocity_context.h"
#include "contexts/thermostat_context.h"

#include <memory>
#include <vector>
#include <string>

class SystemState;
class ForceManager;
struct SimulationConfig;

enum class DumpType : std::uint8_t {
    POSITION,
    VELOCITY,
    FORCE,
    UNKNOWN
};

struct DumpItem {
    std::string name;
    std::string unit;

    [[nodiscard]] bool isEnabled() const {
        return (unit != "off" && unit != "false");
    }

    [[nodiscard]] std::string getEffectiveUnit() const {
        if (unit == "true" || unit == "on") return "atomic_unit";  // Default unit
        if (unit == "none") return "";                             // No unit conversion
        return unit;                                               // Use specified unit
    }

    [[nodiscard]] DumpType getType() const {
        if (name == "positions") return DumpType::POSITION;
        if (name == "velocities") return DumpType::VELOCITY;
        if (name == "forces") return DumpType::FORCE;
        return DumpType::UNKNOWN;
    }
};

class DumpInitializer {
public:
    DumpInitializer(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state,
        const VelocityContext& velocity_context
    );

    // Main interface method
    [[nodiscard]] std::vector<std::shared_ptr<Dump>> createDumps() const;

private:
    [[nodiscard]] std::vector<DumpItem> parseDumpsList() const;

    // Dump creation methods
    [[nodiscard]] std::shared_ptr<Dump> createPositionsDump(const std::string& out_unit) const;
    [[nodiscard]] std::shared_ptr<Dump> createVelocitiesDump(const std::string& out_unit) const;
    [[nodiscard]] std::shared_ptr<Dump> createForcesDump(const std::string& out_unit) const;

    std::shared_ptr<SimulationConfig> m_config;
    std::shared_ptr<SystemState> m_state;

    VelocityContext m_velocity_context;
    ThermostatContext m_thermostat_context;
};
