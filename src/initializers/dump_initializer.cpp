#include "initializers/dump_initializer.h"
#include "core/simulation_config.h"
#include "core/system_state.h"

#include <stdexcept>
#include <memory>

DumpInitializer::DumpInitializer(
    const std::shared_ptr<SimulationConfig>& config,
    const std::shared_ptr<SystemState>& state,
    const VelocityContext& velocity_context
) : m_config(config),
    m_state(state),
    m_velocity_context(velocity_context)
{
}

std::shared_ptr<Dump> DumpInitializer::createPositionsDump(const std::string& out_unit) const {
    return std::make_shared<PositionDump>(
        std::shared_ptr<VecArray>(m_state, &m_state->coord),
        m_config->this_bead,
        m_config->sfreq, /// TODO: Generalize the save frequency option to dumps, observables, etc.
        out_unit
    );
}

std::shared_ptr<Dump> DumpInitializer::createVelocitiesDump(const std::string& out_unit) const {
    return std::make_shared<VelocityDump>(
        m_velocity_context,
        m_config->this_bead,
        m_config->sfreq, /// TODO: Generalize the save frequency option to dumps, observables, etc.
        out_unit
    );
}

std::shared_ptr<Dump> DumpInitializer::createForcesDump(const std::string& out_unit) const {
    return std::make_shared<ForceDump>(
        m_state,
        m_config->this_bead,
        m_config->sfreq, /// TODO: Generalize the save frequency option to dumps, observables, etc.
        out_unit
    );
}

std::vector<DumpItem> DumpInitializer::parseDumpsList() const {
    std::vector<DumpItem> items;
    items.reserve(m_config->dumps_list.size());

    for (const auto& [name, unit] : m_config->dumps_list) {
        //items.emplace_back(DumpItem{ .name = name, .unit = unit });
        items.emplace_back(name, unit);
    }

    return items;
}

std::vector<std::shared_ptr<Dump>> DumpInitializer::createDumps() const {
    const auto items = parseDumpsList();

    std::vector<std::shared_ptr<Dump>> dumps;
    dumps.reserve(items.size());

    for (const auto& item : items) {
        // Skip disabled observables
        if (!item.isEnabled()) {
            continue;
        }

        // Direct method calls based on observable type
        try {
            switch (item.getType()) {
            case DumpType::POSITION:
                dumps.push_back(createPositionsDump(item.getEffectiveUnit()));
                break;
            case DumpType::VELOCITY:
                dumps.push_back(createVelocitiesDump(item.getEffectiveUnit()));
                break;
            case DumpType::FORCE:
                dumps.push_back(createForcesDump(item.getEffectiveUnit()));
                break;
            case DumpType::UNKNOWN:
            default:
                throw std::runtime_error("Unknown dump type: " + item.name);
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to create dump '" + item.name + "': " + e.what());
        }
    }

    return dumps;
}
