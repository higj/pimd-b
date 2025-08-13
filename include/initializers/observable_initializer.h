#pragma once

#include "observable_item.h"
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
class ForceManager;
class Thermostat;
struct SimulationConfig;

class ObservableInitializer
{
public:
    ObservableInitializer(
        const long stride,
        const StringMap& obs_list,
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<ForceManager>& force_mgr,
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

    long m_stride;
    StringMap m_observables_list;
    std::shared_ptr<SystemState> m_state;
    std::shared_ptr<ForceManager> m_force_mgr;

    BeadContext m_bead_context;
    ThermalContext m_thermal_context;
    SpringContext m_spring_context;
    BoxContext m_box_context;
    VelocityContext m_velocity_context;
    ThermostatContext m_thermostat_context;
};
