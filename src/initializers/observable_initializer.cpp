#include "initializers/observable_initializer.h"
#include "core/simulation_config.h"
#include "core/system_state.h"

#include <stdexcept>
#include <memory>

ObservableInitializer::ObservableInitializer(
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
) : m_config(config),
    m_state(state),
    m_exchange_state(exchange_state),
    m_force_mgr(force_mgr),
    m_thermostat(thermostat),
    m_bead_context(bead_context),
    m_thermal_context(thermal_context),
    m_spring_context(spring_context),
    m_box_context(box_context),
    m_velocity_context(velocity_context),
    m_thermostat_context(thermostat_context)
{
}

std::shared_ptr<Observable> ObservableInitializer::createEnergyObservable(const std::string& out_unit) const {

    return std::make_shared<EnergyObservable>(
        m_exchange_state,
        std::shared_ptr<const dVec>(m_state, &m_state->coord),
        std::shared_ptr<const dVec>(m_state, &m_state->prev_coord),
        m_force_mgr,
        m_bead_context,
        m_thermal_context,
        m_spring_context,
        m_box_context,
        m_config->sfreq,
        out_unit
    );
}

std::shared_ptr<Observable> ObservableInitializer::createClassicalObservable(const std::string& out_unit) const {

    return std::make_shared<ClassicalObservable>(
        std::shared_ptr<const dVec>(m_state, &m_state->coord),
        std::shared_ptr<const dVec>(m_state, &m_state->prev_coord),
        m_exchange_state,
        m_velocity_context,
        m_thermostat_context,
        m_bead_context,
        m_spring_context,
        m_box_context,
        m_config->sfreq,
        out_unit
    );
}

std::shared_ptr<Observable> ObservableInitializer::createBosonicObservable(const std::string& out_unit) const {

    // Validate configuration
    if (!m_config->bosonic) {
        throw std::runtime_error("Bosonic observables require bosonic simulation mode");
    }

    return std::make_shared<BosonicObservable>(
        m_exchange_state,
        m_config->sfreq,
        out_unit
    );
}

std::shared_ptr<Observable> ObservableInitializer::createGSFObservable(const std::string& out_unit) const {
    return std::make_shared<GSFActionObservable>(
        std::shared_ptr<const dVec>(m_state, &m_state->coord),
        m_force_mgr,
        m_bead_context,
        m_thermal_context,
        m_spring_context,
        m_config->sfreq,
        out_unit
    );
}

std::vector<ObservableItem> ObservableInitializer::parseObservablesList() const {
    std::vector<ObservableItem> items;
    items.reserve(m_config->observables_list.size());

    for (const auto& [name, unit] : m_config->observables_list) {
        items.emplace_back(ObservableItem{.name = name, .unit = unit });
    }

    return items;
}

std::vector<std::shared_ptr<Observable>> ObservableInitializer::createObservables() const {
    const auto items = parseObservablesList();

    std::vector<std::shared_ptr<Observable>> observables;
    observables.reserve(items.size());

    for (const auto& item : items) {
        // Skip disabled observables
        if (!item.isEnabled()) {
            continue;
        }

        // Direct method calls based on observable name
        try {
            if (item.name == "energy") {
                observables.push_back(createEnergyObservable(item.getEffectiveUnit()));
            } else if (item.name == "classical") {
                observables.push_back(createClassicalObservable(item.getEffectiveUnit()));
            } else if (item.name == "bosonic") {
                observables.push_back(createBosonicObservable(item.getEffectiveUnit()));
            } else if (item.name == "gsf") {
                observables.push_back(createGSFObservable(item.getEffectiveUnit()));
            } else {
                throw std::runtime_error("Unknown observable type: " + item.name);
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to create observable '" + item.name + "': " + e.what());
        }
    }

    return observables;
}