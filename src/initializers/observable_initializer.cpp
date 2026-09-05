#include "initializers/observable_initializer.h"
#include "observables.h"
#include "core/simulation_config.h"
#include "core/system_state.h"

#include <stdexcept>
#include <memory>
#include <string_view>

namespace {
    // Output destinations are deliberately defined in code, not in the input file.
    // Leave an entry empty to write that observable to Output::MAIN_FILENAME.
    constexpr std::string_view ENERGY_OUTPUT_FILENAME;
    constexpr std::string_view CLASSICAL_OUTPUT_FILENAME;
    constexpr std::string_view BOSONIC_OUTPUT_FILENAME;
    constexpr std::string_view GSF_OUTPUT_FILENAME;
    constexpr std::string_view CENTER_OF_MASS_OUTPUT_FILENAME;
}

ObservableInitializer::ObservableInitializer(
    const long stride,
    const StringMap& obs_list,
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<ForceManager>& force_mgr,
    const ActionContext& action_context,
    const BeadContext& bead_context,
    const ThermalContext& thermal_context,
    const SpringContext& spring_context,
    const BoxContext& box_context,
    const VelocityContext& velocity_context,
    const ThermostatContext& thermostat_context
) : m_stride(stride),
    m_observables_list(std::move(obs_list)),
    m_state(state),
    m_force_mgr(force_mgr),
    m_action_context(action_context),
    m_bead_context(bead_context),
    m_thermal_context(thermal_context),
    m_spring_context(spring_context),
    m_box_context(box_context),
    m_velocity_context(velocity_context),
    m_thermostat_context(thermostat_context)
{
}

std::shared_ptr<Observable> ObservableInitializer::createPrimitiveKineticEnergyObservable(
    const std::string& out_unit) const {
    auto obs = std::make_shared<PrimitiveKineticEnergyObservable>(
        std::shared_ptr<const VecArray>(m_state, &m_state->coord),
        std::shared_ptr<const VecArray>(m_state, &m_state->prev_coord),
        m_bead_context, m_thermal_context, m_spring_context, m_box_context,
        m_stride, out_unit);
    obs->setOutputFilename(ENERGY_OUTPUT_FILENAME);
    return obs;
}

std::shared_ptr<Observable> ObservableInitializer::createPotentialEnergyObservable(
    const std::string& out_unit) const {
    auto obs = std::make_shared<PotentialEnergyObservable>(
        std::shared_ptr<const VecArray>(m_state, &m_state->coord),
        m_force_mgr, m_box_context, m_bead_context,
        m_stride, out_unit);
    obs->setOutputFilename(ENERGY_OUTPUT_FILENAME);
    return obs;
}

std::shared_ptr<Observable> ObservableInitializer::createVirialEnergyObservable(
    const std::string& out_unit) const {
    auto obs = std::make_shared<VirialEnergyObservable>(
        std::shared_ptr<const VecArray>(m_state, &m_state->coord),
        std::shared_ptr<const VecArray>(m_state, &m_state->physical_forces),
        m_bead_context, m_stride, out_unit);
    obs->setOutputFilename(ENERGY_OUTPUT_FILENAME);
    return obs;
}

std::shared_ptr<Observable> ObservableInitializer::createExtPotObservable(
    const std::string& out_unit) const {
    auto obs = std::make_shared<ExtPotObservable>(
        std::shared_ptr<const VecArray>(m_state, &m_state->coord),
        m_force_mgr, m_bead_context, m_stride, out_unit);
    obs->setOutputFilename(ENERGY_OUTPUT_FILENAME);
    return obs;
}

std::shared_ptr<Observable> ObservableInitializer::createIntPotObservable(
    const std::string& out_unit) const {
    auto obs = std::make_shared<IntPotObservable>(
        std::shared_ptr<const VecArray>(m_state, &m_state->coord),
        m_force_mgr, m_box_context, m_bead_context, m_stride, out_unit);
    obs->setOutputFilename(ENERGY_OUTPUT_FILENAME);
    return obs;
}

std::shared_ptr<Observable> ObservableInitializer::createClassicalKineticEnergyObservable(const std::string& out_unit) const
{
    auto observable = std::make_shared<ClassicalKineticEnergyObservable>(
        m_velocity_context,
        m_bead_context,
        m_stride,
        out_unit
    );
    observable->setOutputFilename(CLASSICAL_OUTPUT_FILENAME);
    return observable;
}

std::shared_ptr<Observable> ObservableInitializer::createClassicalSpringEnergyObservable(const std::string& out_unit) const
{
    auto observable = std::make_shared<ClassicalSpringEnergyObservable>(
        std::shared_ptr<const VecArray>(m_state, &m_state->coord),
        std::shared_ptr<const VecArray>(m_state, &m_state->prev_coord),
        m_spring_context,
        m_box_context,
        m_stride,
        out_unit
    );
    observable->setOutputFilename(CLASSICAL_OUTPUT_FILENAME);
    return observable;
}

std::shared_ptr<Observable> ObservableInitializer::createTemperatureObservable(const std::string& out_unit) const
{
    auto observable = std::make_shared<TemperatureObservable>(
        m_velocity_context,
        m_bead_context,
        m_stride,
        out_unit
    );
    observable->setOutputFilename(CLASSICAL_OUTPUT_FILENAME);
    return observable;
}

std::shared_ptr<Observable> ObservableInitializer::createNoseHooverEnergyObservable(const std::string& out_unit) const {
    auto observable = std::make_shared<NoseHooverEnergyObservable>(
        m_thermostat_context,
        m_stride,
        out_unit
    );
    observable->setOutputFilename(CLASSICAL_OUTPUT_FILENAME);
    return observable;
}

std::shared_ptr<Observable> ObservableInitializer::createBosonicProbDistObservable(const std::string& out_unit) const
{
    if (!m_state->isBosonic())
    {
        throw std::runtime_error("Bosonic observables require bosonic simulation mode");
    }

    auto observable = std::make_shared<BosonicProbDistObservable>(
        m_stride,
        out_unit
    );
    observable->setOutputFilename(BOSONIC_OUTPUT_FILENAME);
    return observable;
}

std::shared_ptr<Observable> ObservableInitializer::createBosonicProbAllObservable(const std::string& out_unit) const
{
    if (!m_state->isBosonic())
    {
        throw std::runtime_error("Bosonic observables require bosonic simulation mode");
    }
    auto observable = std::make_shared<BosonicProbAllObservable>(
        m_stride,
        out_unit
    );
    observable->setOutputFilename(BOSONIC_OUTPUT_FILENAME);
    return observable;
}

std::shared_ptr<Observable> ObservableInitializer::createBosonicSignObservable(const std::string& out_unit) const
{
    if (!m_state->isBosonic())
    {
        throw std::runtime_error("Bosonic observables require bosonic simulation mode");
    }
    auto observable = std::make_shared<BosonicSignObservable>(
        m_stride,
        out_unit
    );
    observable->setOutputFilename(BOSONIC_OUTPUT_FILENAME);
    return observable;
}

std::shared_ptr<Observable> ObservableInitializer::createSuzukiChinWeightObservable(
    const std::string& out_unit) const {
    auto obs = std::make_shared<SuzukiChinWeightObservable>(
        std::shared_ptr<const VecArray>(m_state, &m_state->coord),
        m_force_mgr, m_bead_context, m_box_context, m_thermal_context, m_spring_context,
        m_action_context.gsf_alpha, m_stride, out_unit);
    obs->setOutputFilename(GSF_OUTPUT_FILENAME);
    return obs;
}

std::shared_ptr<Observable> ObservableInitializer::createSuzukiChinPotEnergyObservable(
    const std::string& out_unit) const {
    auto obs = std::make_shared<SuzukiChinPotEnergyObservable>(
        std::shared_ptr<const VecArray>(m_state, &m_state->coord),
        m_force_mgr, m_bead_context, m_box_context, m_stride, out_unit);
    obs->setOutputFilename(GSF_OUTPUT_FILENAME);
    return obs;
}

std::shared_ptr<Observable> ObservableInitializer::createSuzukiChinEvenPotEnergyObservable(
    const std::string& out_unit) const {
    auto obs = std::make_shared<SuzukiChinEvenPotEnergyObservable>(
        std::shared_ptr<const VecArray>(m_state, &m_state->coord),
        m_force_mgr, m_bead_context, m_box_context, m_stride, out_unit);
    obs->setOutputFilename(GSF_OUTPUT_FILENAME);
    return obs;
}

std::shared_ptr<Observable> ObservableInitializer::createSuzukiChinKineticEnergyObservable(
    const std::string& out_unit) const {
    auto obs = std::make_shared<SuzukiChinKineticEnergyObservable>(
        std::shared_ptr<const VecArray>(m_state, &m_state->coord),
        m_force_mgr, m_bead_context, m_box_context, m_spring_context,
        m_action_context.gsf_alpha, m_stride, out_unit);
    obs->setOutputFilename(GSF_OUTPUT_FILENAME);
    return obs;
}

std::shared_ptr<Observable> ObservableInitializer::createSuzukiChinVirialEnergyObservable(
    const std::string& out_unit) const {
    auto obs = std::make_shared<SuzukiChinVirialEnergyObservable>(
        std::shared_ptr<const VecArray>(m_state, &m_state->coord),
        m_force_mgr, m_bead_context, m_box_context, m_stride, out_unit);
    obs->setOutputFilename(GSF_OUTPUT_FILENAME);
    return obs;
}

std::shared_ptr<Observable> ObservableInitializer::createCenterOfMassObservable(const std::string& out_unit) const
{
    auto observable = std::make_shared<CenterOfMassObservable>(
        std::shared_ptr<const VecArray>(m_state, &m_state->coord),
        m_bead_context,
        m_stride,
        out_unit
    );
    observable->setOutputFilename(CENTER_OF_MASS_OUTPUT_FILENAME);
    return observable;
}

std::vector<ObservableItem> ObservableInitializer::parseObservablesList() const
{
    std::vector<ObservableItem> items;
    items.reserve(m_observables_list.size());

    for (const auto& [name, unit] : m_observables_list)
    {
        //items.emplace_back(ObservableItem{.name = name, .unit = unit});
        items.emplace_back(name, unit);
    }

    return items;
}

std::vector<std::shared_ptr<Observable>> ObservableInitializer::createObservables() const
{
    const auto items = parseObservablesList();

    std::vector<std::shared_ptr<Observable>> observables;
    observables.reserve(items.size());

    for (const auto& item : items)
    {
        // Skip disabled observables
        if (!item.isEnabled())
        {
            continue;
        }

        // Direct method calls based on observable type
        try
        {
            switch (item.getType())
            {
            case ObservableType::PRIM_KE:
                observables.push_back(createPrimitiveKineticEnergyObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::VIRIAL_KE:
                observables.push_back(createVirialEnergyObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::POT_ENERGY:
                observables.push_back(createPotentialEnergyObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::EXT_POT:
                observables.push_back(createExtPotObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::INT_POT:
                 observables.push_back(createIntPotObservable(item.getEffectiveUnit()));
                 break;
            case ObservableType::CL_KINETIC_ENERGY:
                observables.push_back(createClassicalKineticEnergyObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::CL_SPRING_ENERGY:
                observables.push_back(createClassicalSpringEnergyObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::QUANTUM_TEMPERATURE:
                observables.push_back(createTemperatureObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::NOSE_HOOVER_ENERGY:
                observables.push_back(createNoseHooverEnergyObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::BOSONIC_PROB_DIST:
                observables.push_back(createBosonicProbDistObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::BOSONIC_PROB_ALL:
                observables.push_back(createBosonicProbAllObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::BOSONIC_SIGN:
                observables.push_back(createBosonicSignObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::SC_WEIGHT:
                observables.push_back(createSuzukiChinWeightObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::SC_POT_ENERGY:
                observables.push_back(createSuzukiChinPotEnergyObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::SC_EVEN_POT_ENERGY:
                observables.push_back(createSuzukiChinEvenPotEnergyObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::SC_KINETIC_ENERGY:
                observables.push_back(createSuzukiChinKineticEnergyObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::SC_VIRIAL_ENERGY:
                observables.push_back(createSuzukiChinVirialEnergyObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::CENTER_OF_MASS:
                observables.push_back(createCenterOfMassObservable(item.getEffectiveUnit()));
                break;
            case ObservableType::UNKNOWN:
            default:
                throw std::runtime_error("Unknown observable type: " + item.name);
            }
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error("Failed to create observable '" + item.name + "': " + e.what());
        }
    }

    return observables;
}
