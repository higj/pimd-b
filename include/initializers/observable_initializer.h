#pragma once

#include "observables.h"
#include "contexts/bead_context.h"
#include "contexts/thermal_context.h"
#include "contexts/spring_context.h"
#include "contexts/box_context.h"
#include "contexts/velocity_context.h"
#include "contexts/thermostat_context.h"

#include <memory>
#include <vector>
#include <string>

class SystemState;
class ExchangeState;
class ForceManager;
class Thermostat;
struct SimulationConfig;

enum class ObservableType : std::uint8_t {
    ENERGY,
    CLASSICAL,
    BOSONIC,
    GSF,
    UNKNOWN
};

struct ObservableItem
{
    std::string name;
    std::string unit;

    [[nodiscard]] bool isEnabled() const
    {
        return unit != "off";
    }

    // If the unit is "none", return an empty string, because no unit conversion will take place
    [[nodiscard]] std::string getEffectiveUnit() const
    {
        return (unit != "none") ? unit : "";
    }

    [[nodiscard]] ObservableType getType() const
    {
        if (name == "energy") return ObservableType::ENERGY;
        if (name == "classical") return ObservableType::CLASSICAL;
        if (name == "bosonic") return ObservableType::BOSONIC;
        if (name == "gsf") return ObservableType::GSF;
        return ObservableType::UNKNOWN;
    }
};

class ObservableInitializer
{
public:
    ObservableInitializer(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<ExchangeState>& exchange_state,
        const std::shared_ptr<ForceManager>& force_mgr,
        const std::shared_ptr<Thermostat>& thermostat,
        const BeadContext& bead_context,
        const ThermalContext& thermal_context,
        const SpringContext& spring_context,
        const BoxContext& box_context,
        const VelocityContext& velocity_context,
        const ThermostatContext& thermostat_context
    );

    // Main interface method
    [[nodiscard]] std::vector<std::shared_ptr<Observable>> createObservables() const;

private:
    [[nodiscard]] std::vector<ObservableItem> parseObservablesList() const;

    // Observable creation methods
    [[nodiscard]] std::shared_ptr<Observable> createEnergyObservable(const std::string& out_unit) const;
    [[nodiscard]] std::shared_ptr<Observable> createClassicalObservable(const std::string& out_unit) const;
    [[nodiscard]] std::shared_ptr<Observable> createBosonicObservable(const std::string& out_unit) const;
    [[nodiscard]] std::shared_ptr<Observable> createGSFObservable(const std::string& out_unit) const;

    std::shared_ptr<SimulationConfig> m_config;
    std::shared_ptr<SystemState> m_state;
    std::shared_ptr<ExchangeState> m_exchange_state;
    std::shared_ptr<ForceManager> m_force_mgr;
    std::shared_ptr<Thermostat> m_thermostat;

    BeadContext m_bead_context;
    ThermalContext m_thermal_context;
    SpringContext m_spring_context;
    BoxContext m_box_context;
    VelocityContext m_velocity_context;
    ThermostatContext m_thermostat_context;
};
