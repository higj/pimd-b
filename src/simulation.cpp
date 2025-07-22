#include "simulation.h"
#include "params.h"
#include "core/system_state.h"
#include "core/force_manager.h"
#include "core/random_generators.h"
#include "core/statistics_manager.h"
#include "momentum_initializers.h"
#include "position_initializers.h"
#include "initializers/observable_initializer.h"
#include "initializers/dump_initializer.h"
#include "propagators.h"
#include "thermostats.h"
#include "observables_logger.h"
#include "output_paths.h"
#include "simulation_report.h"

#include <filesystem>
#include <array>
#include <cassert>

Simulation::Simulation(int rank, int nproc, const std::string& config_filename)
{
    // Parse the simulation parameters from the configuration (input) file
    const Params params(config_filename, rank);
    const std::shared_ptr<SimulationConfig> config = params.load();

    m_steps = config->steps;
    m_threshold = config->threshold;

    // Initialize basic contexts
    m_thermal_ctx = ThermalContext{
        .beta = config->beta,
        .thermo_beta = config->thermo_beta
    };

    m_spring_ctx = SpringContext{
        .omega_p = config->omega_p,
        .spring_constant = config->spring_constant,
        .beta_half_k = config->beta_half_k
    };

    m_box_ctx = BoxContext{
        .box_size = config->box_size,
        .pbc = config->pbc
    };

    m_bead_ctx = BeadContext{
        .nbeads = config->nbeads,
        .natoms = config->natoms,
        .this_bead = config->this_bead
    };

    // CR: What is meaning of the variable prefix m_?
    // JH: Hungarian notation for member of a class.
    //     Arguably better than prefixing with just an underscore, or using "this->"
    //     to distinguish from local variables.
    m_rng = std::make_shared<RandomGenerators>(config->seed + rank);
    m_state = std::make_shared<SystemState>(rank, nproc, config->natoms, config->nbeads, config->fixcom);

    initializePositions(config);
    initializeMomenta(config);

    initializeQuantumStatistics(config->bosonic);

    // Initialize other resources
    m_force_mgr = std::make_shared<ForceManager>(*config);    

    m_normal_modes = initializeNormalModes(config, m_state);

    m_nm_ctx = NormalModesContext{
        .normal_modes = m_normal_modes,
        .couple_to_nm = std::get<bool>(config->thermostat_params["nmthermostat"])
    };

    // CR: It's OK imo to pass all of config to the local function (of the class Simulation) that initializes the propagator.
    // CR: Would improve readability of this function.
    // CR: The important bit is that Config & Simulation are not exposed to *other* classes.
    m_propagator = initializePropagator(
        config,
        m_state,
        m_normal_modes,
        m_force_mgr
    );

    m_thermostat = initializeThermostat(
        config,
        m_state,
        m_rng
    );

    m_thermostat_ctx = ThermostatContext{
        .thermostat = m_thermostat,
        .thermostat_type = config->thermostat_type
    };

    m_velocity_ctx = VelocityContext{
        .momenta = std::shared_ptr<dVec>(m_state, &m_state->momenta),
        .mass = config->mass
    };

    /// TODO: Should this also be done in dumps and observables init for safety?
    std::filesystem::create_directory(Output::FOLDER_NAME);

    // CR: I'd write this->m_observables for emphasis
    // CR: I think it's OK for these classes to accept config too, because they have
    // CR: no other logic except bridging the gap between config and the simulation.
    // CR: They don't continue to live after initialization
    m_observables = ObservableInitializer(
        config,
        m_state,
        m_force_mgr,
        m_thermostat,
        m_bead_ctx,
        m_thermal_ctx,
        m_spring_ctx,
        m_box_ctx,
        m_velocity_ctx,
        m_thermostat_ctx
    ).createObservables();

    // Initialize the output file for the observables
    m_obs_logger = std::make_unique<ObservablesLogger>(
        Output::MAIN_FILENAME, m_bead_ctx.this_bead, config->sfreq, m_observables
    );

    m_dumps = initializeDumps(config, m_velocity_ctx);
    m_report = std::make_unique<SimulationReport>(*config, m_steps);
}

Simulation::~Simulation() = default;

void Simulation::initializeQuantumStatistics(bool is_bosonic) const
{
    if (!m_state)
    {
        throw std::runtime_error("System state is not initialized before quantum statistics initialization.");
    }

    StatisticsManager::getInstance().initializeBosonic(
        is_bosonic,
        m_bead_ctx,
        m_thermal_ctx,
        m_spring_ctx,
        m_box_ctx,
        m_state
    );
}

std::shared_ptr<NormalModes> Simulation::initializeNormalModes(
    const std::shared_ptr<SimulationConfig>& config,
    const std::shared_ptr<SystemState>& state)
{
    /// TODO: Perhaps we should initialize NormalModes regardless of the propagator type and/or coupling?
    if (const bool couple_to_nm = std::get<bool>(config->thermostat_params["nmthermostat"]);
        config->propagator_type == "normal_modes" || couple_to_nm)
    {
        return std::make_shared<NormalModes>(
            std::shared_ptr<dVec>(state, &state->coord),
            std::shared_ptr<dVec>(state, &state->momenta),
            config->natoms,
            config->nbeads,
            config->this_bead
        );
    }
    return {nullptr};
}

std::shared_ptr<Propagator> Simulation::initializePropagator(
    const std::shared_ptr<SimulationConfig>& config,
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<NormalModes>& normal_modes,
    const std::shared_ptr<ForceManager>& force_mgr
)
{
    if (config->propagator_type == "cartesian")
    {
        return std::make_shared<VelocityVerletPropagator>(
            state,
            force_mgr,
            m_box_ctx,
            m_spring_ctx,
            config->mass,
            config->dt
        );
    }

    if (config->propagator_type == "normal_modes")
    {
        return std::make_shared<NormalModesPropagator>(
            state,
            force_mgr,
            m_box_ctx,
            m_spring_ctx,
            config->mass,
            config->dt,
            normal_modes
        );
    }

    return {nullptr};
}

std::shared_ptr<Thermostat> Simulation::initializeThermostat(
    const std::shared_ptr<SimulationConfig>& config,
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<RandomGenerators>& rng)
{
    if (config->thermostat_type == "langevin")
    {
        double gamma = std::get<double>(config->thermostat_params["gamma"]);

        return std::make_shared<LangevinThermostat>(
            m_thermal_ctx,
            m_nm_ctx,
            state,
            rng,
            gamma,
            config->dt,
            config->mass
        );
    }

    int nchains = std::get<int>(config->thermostat_params["nchains"]);

    if (config->thermostat_type == "nose_hoover")
    {
        return std::make_shared<NoseHooverThermostat>(m_thermal_ctx, m_nm_ctx, state, nchains, config->dt, config->mass);
    }

    if (config->thermostat_type == "nose_hoover_np")
    {
        return std::make_shared<NoseHooverNpThermostat>(m_thermal_ctx, m_nm_ctx, state, nchains, config->dt, config->mass);
    }

    if (config->thermostat_type == "nose_hoover_np_dim")
    {
        return std::make_shared<NoseHooverNpDimThermostat>(m_thermal_ctx, m_nm_ctx, state, nchains, config->dt,
                                                           config->mass);
    }

    if (config->thermostat_type == "none")
    {
        return std::make_shared<Thermostat>(m_thermal_ctx, m_nm_ctx, state);
    }

    return {nullptr};
}

void Simulation::initializePositions(const std::shared_ptr<SimulationConfig>& config)
{
    if (!m_state || !m_rng) {
        throw std::runtime_error("State or RNG is not initialized before position initialization.");
    }

    std::unique_ptr<PositionInitializer> initializer;

    if (config->init_pos_type == "random")
    {
        initializer = std::make_unique<RandomPositionInitializer>(
            m_rng,
            std::shared_ptr<dVec>(m_state, &m_state->coord),
            m_box_ctx
        );
    }
    else if (config->init_pos_type == "grid")
    {
        initializer = std::make_unique<GridPositionInitializer>(
            std::shared_ptr<dVec>(m_state, &m_state->coord),
            m_box_ctx
        );
    }
    else if (config->init_pos_type == "xyz")
    {
        initializer = std::make_unique<XyzPositionInitializer>(
            config->init_pos_filename,
            config->this_bead + config->init_pos_index_offset,
            std::shared_ptr<dVec>(m_state, &m_state->coord),
            m_box_ctx
        );
    }
    else
    {
        throw std::invalid_argument("Unknown position initialization method: " + config->init_pos_type);
    }

    initializer->initialize();

    // Communicate the new coordinates to the neighboring processes
    m_state->updateNeighboringCoordinates();
}

void Simulation::initializeMomenta(const std::shared_ptr<SimulationConfig>& config)
{
    if (!m_state || !m_rng)
    {
        throw std::runtime_error("State or RNG is not initialized before momentum initialization.");
    }

    std::unique_ptr<MomentumInitializer> initializer;

    if (config->init_vel_type == "random")
    {
        initializer = std::make_unique<MaxwellBoltzmannMomentumInitializer>(
            m_rng,
            m_state,
            config->mass,
            config->thermo_beta
        );
    }
    else if (config->init_vel_type == "manual")
    {
        initializer = std::make_unique<ManualMomentumInitializer>(
            config->init_vel_filename,
            config->init_vel_index_offset,
            m_state,
            config->mass
        );
    }
    else
    {
        throw std::invalid_argument("Unknown momentum initialization method: " + config->init_vel_type);
    }

    initializer->initialize();
}

void Simulation::setStep(long step)
{
    m_is_thermalization_phase = (step < m_threshold);
}

std::vector<std::shared_ptr<Dump>> Simulation::initializeDumps(
    const std::shared_ptr<SimulationConfig>& config, 
    const VelocityContext& vel_ctx
) const
{
    auto dumps = DumpInitializer(config, m_state, vel_ctx).createDumps();

    for (const auto& dump : dumps) {
        dump->initialize();
    }

    return dumps;
}

double Simulation::getWallTime()
{
    return MPI_Wtime();
}

void Simulation::zeroMomentumIfRequired() const
{
    if (m_state->isCenterOfMassFixed()) {
        m_state->zeroMomentum();
    }
}

void Simulation::resetObservables() const {
    for (const auto& observable : m_observables) {
        observable->resetValues();
    }
}

void Simulation::dumpStepInfo(long step) const
{
    for (const auto& dump : m_dumps) {
        dump->output(step);
    }
}

void Simulation::performMolecularDynamicsStep() const {
    m_thermostat->step();
    zeroMomentumIfRequired();
    m_propagator->step();
    m_thermostat->step();
    zeroMomentumIfRequired();
}

void Simulation::calculateObservables() const {
    for (const auto& observable : m_observables) {
        observable->calculate();
    }
}

void Simulation::calculateAndLogObservables(long step) const {
    calculateObservables();
    m_obs_logger->log(step);
}

void Simulation::finalizeSimulation(double start_time) const {
    const double wall_time = getWallTime() - start_time;

    printStatus(
        std::format(
            "Simulation finished running successfully (Runtime = {:.3} sec)", 
            wall_time
        ),
        m_bead_ctx.this_bead
    );

    m_report->writeReport(wall_time);
}

void Simulation::executeStep(long step) const
{
    resetObservables();
    dumpStepInfo(step);
    performMolecularDynamicsStep();

    if (!m_is_thermalization_phase) {
        calculateAndLogObservables(step);
    }
}

void Simulation::run()
{
    printStatus("Running the simulation", m_bead_ctx.this_bead);
    const double start_time = getWallTime();
    
    // Main loop performing molecular dynamics steps
    for (long step = 0; step <= m_steps; ++step)
    {
        setStep(step);
        executeStep(step);
    }

    finalizeSimulation(start_time);
}