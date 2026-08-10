#include "simulation.h"
#include "mpi_utils.h"
#include "core/simulation_config.h"
#include "core/system_state.h"
#include "core/force_manager.h"
#include "core/random_generators.h"
#include "core/statistics_manager.h"
#include "momentum_initializers.h"
#include "position_initializers.h"
#include "initializers/thermostat_initializer.h"
#include "initializers/observable_initializer.h"
#include "initializers/dump_initializer.h"
#include "initializers/rpmd_frame_selector.h"
#include "propagators.h"
#include "observables_logger.h"
#include "output_paths.h"
#include "simulation_report.h"
#include "thermostats/thermostat.h"

#include <filesystem>
#include <array>
#include <cassert>

// CR: What is meaning of the variable prefix m_?
// JH: Hungarian notation for member of a class.
//     Arguably better than prefixing with just an underscore, or using "this->"
//     to distinguish from local variables.
Simulation::Simulation(int rank, int nproc, const std::shared_ptr<SimulationConfig>& config)
    : m_steps(config->steps),
      m_threshold(config->threshold),
      m_sfreq(config->sfreq),
      m_is_rpmd(config->rpmd_config.enabled),
      m_is_thermalization_phase(true),
      m_state(std::make_shared<SystemState>(rank, nproc, config->natoms, config->nbeads, config->fixcom, config->bosonic)),
      m_rng(std::make_shared<RandomGenerators>(config->seed + rank))
{
    initializeConfigurationDependentContexts(config);
    initializePositions(config);
    initializeMomenta(config);
    initializeQuantumStatistics(config->bosonic);

    // Initialize RPMD frame selector if RPMD is enabled. This will determine which frame to use for each independent run
    if (m_is_rpmd)
    {
        m_rpmd_frame_selector = std::make_unique<RpmdFrameSelector>(config);
    }

    // createPotential() is called exactly once per config; it moves ownership
    // into ForceManager and will throw if accidentally called again.
    m_force_mgr = std::make_shared<ForceManager>(
        config->ext_potential_cfg.createPotential(),
        config->int_potential_cfg.createPotential(),
        config->int_potential_cfg.cutoff(),
        m_box_ctx
    );
    //initializeForces(m_state, m_force_mgr, m_spring_ctx, m_box_ctx); <--- TODO: Uncomment this later

    m_normal_modes = initializeNormalModes(config, m_state);

    m_nm_ctx = NormalModesContext{
        .normal_modes = m_normal_modes,
        // TODO: Here we assume that "nmthermostat" exists. Currently it does, but what if thermostat is "none"?
        .couple_to_nm = config->thermostat.couple_to_nm
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

    /// TODO: Note that config is passed here to a function outside of the Simulation class. Is this okay?
    m_thermostat = ThermostatInitializer(
        config,
        m_thermal_ctx,
        m_nm_ctx
    ).createThermostat(m_state, m_rng);

    m_thermostat_ctx = ThermostatContext{
        .thermostat = m_thermostat,
        .thermostat_type = config->thermostat.type
    };

    m_velocity_ctx = VelocityContext{
        .momenta = std::shared_ptr<VecArray>(m_state, &m_state->momenta),
        .mass = config->mass
    };

    // CR: I'd write this->m_observables for emphasis
    // CR: I think it's OK for these classes to accept config too, because they have
    // CR: no other logic except bridging the gap between config and the simulation.
    // CR: They don't continue to live after initialization
    m_observables = ObservableInitializer(
        config->sfreq,
        config->observables_list,
        m_state,
        m_force_mgr,
        m_bead_ctx,
        m_thermal_ctx,
        m_spring_ctx,
        m_box_ctx,
        m_velocity_ctx,
        m_thermostat_ctx
    ).createObservables();

    /*
    // Initialize the output file for the observables
    m_obs_logger = std::make_unique<ObservablesLogger>(
        Output::MAIN_FILENAME, m_bead_ctx.this_bead, config->sfreq, m_observables
    );
    */

    /// TODO: Deal with dumps and report when RPMD is enabled. Currently, we only support dumps and report for a single run.
    //m_dumps = initializeDumps(config, m_velocity_ctx);
    m_dumps = DumpInitializer(config, m_state, m_velocity_ctx).createDumps();
    m_report = std::make_unique<SimulationReport>(*config, m_steps);
}

Simulation::~Simulation() = default;

void Simulation::initializeConfigurationDependentContexts(
    const std::shared_ptr<SimulationConfig>& config
)
{
    m_rpmd_context = RpmdContext{
        .enabled = config->rpmd_config.enabled,
        .num_runs = config->rpmd_config.num_runs
    };

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
}

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

/// TODO: Refactoring this into a separate class is probably an overkill (NM are either used or not,
///       and there are not too many options here). Perhaps create a separate namespace with a
///       single function that will initialize the normal modes?
std::shared_ptr<NormalModes> Simulation::initializeNormalModes(
    const std::shared_ptr<SimulationConfig>& config,
    const std::shared_ptr<SystemState>& state)
{
    /// TODO: Perhaps we should initialize NormalModes regardless of the propagator type and/or coupling?
    if (const bool couple_to_nm = config->thermostat.couple_to_nm;
        config->propagator_type == "normal_modes" || couple_to_nm)
    {
        return std::make_shared<NormalModes>(
            std::shared_ptr<VecArray>(state, &state->coord),
            std::shared_ptr<VecArray>(state, &state->momenta),
            config->natoms,
            config->nbeads,
            config->this_bead
        );
    }

    return {};
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

    return {};
}

void Simulation::bindPositionInitFactory(
    const std::shared_ptr<SimulationConfig>& config
) {
    const std::string init_type = config->init_pos_type;
    const std::string filename = config->init_pos_filename;
    const int first_idx = config->this_bead + config->init_pos_index_offset;
    const std::string unit = config->init_pos_unit;
    //const XyzFrameSelectionMode frame_mode = frame_selection_mode.value_or(config->init_pos_frame_mode);
    const XyzFrameSelectionMode frame_mode = config->rpmd_config.enabled ? XyzFrameSelectionMode::Index : config->init_pos_frame_mode;

    m_position_init_factory = [this, init_type, filename, first_idx, unit, frame_mode, config]
        (const std::optional<long> override_frame)
    {
        if (!m_state || !m_rng) {
            throw std::runtime_error("State or RNG is not initialized before position initialization.");
        }

        std::unique_ptr<PositionInitializer> initializer;

        if (init_type == "random") {
            initializer = std::make_unique<RandomPositionInitializer>(
                m_rng,
                std::shared_ptr<VecArray>(m_state, &m_state->coord),
                m_box_ctx
            );
        } else if (init_type == "grid") {
            initializer = std::make_unique<GridPositionInitializer>(
                std::shared_ptr<VecArray>(m_state, &m_state->coord),
                m_box_ctx
            );
        } else if (init_type == "xyz") {
            //const long effective_frame = override_frame.has_value() ? override_frame.value() : config->init_pos_frame;
            //const long effective_frame = override_frame.value_or(config->init_pos_frame);
            long effective_frame;
            if (override_frame.has_value()) {
                effective_frame = override_frame.value();
            } else {
                if (config->rpmd_config.enabled) {
                    // If RPMD is enabled, set the effective frame to be the first frame of the XYZ file,
                    // and let the RPMD frame selector handle the frame selection for each run later.
                    effective_frame = 0;
                } else {
                    effective_frame = config->init_pos_frame;
                }
            }

            initializer = std::make_unique<XyzPositionInitializer>(
                filename,
                first_idx,
                unit,
                effective_frame,
                frame_mode,
                std::shared_ptr<VecArray>(m_state, &m_state->coord),
                m_box_ctx
            );
        } else {
            throw std::invalid_argument("Unknown position initialization method: " + init_type);
        }

        MpiUtils::runCollectively(
            "Position initialization",
            [&] {
                initializer->initialize();
            }
        );


        // Communicate the new coordinates to the neighboring processes
        m_state->updateNeighboringCoordinates();
    };
}

void Simulation::bindMomentumInitFactory(
    const std::shared_ptr<SimulationConfig>& config
) {
    const std::string init_type = config->init_vel_type;
    const std::string filename = config->init_vel_filename;
    const int first_idx = config->this_bead + config->init_vel_index_offset;
    const std::string unit = config->init_vel_unit;
    const XyzFrameSelectionMode frame_mode = config->rpmd_config.enabled ? XyzFrameSelectionMode::Index : config->init_vel_frame_mode;
    const double mass = config->mass;
    const double thermo_beta = config->thermo_beta;

    m_momentum_init_factory = [this, init_type, filename, first_idx, unit, frame_mode, mass, thermo_beta, config]
        (const std::optional<long> override_frame)
    {
        if (!m_state || !m_rng) {
            throw std::runtime_error("State or RNG is not initialized before momentum initialization.");
        }

        std::unique_ptr<MomentumInitializer> initializer;

        if (init_type == "random") {
            initializer = std::make_unique<MaxwellBoltzmannMomentumInitializer>(
                m_rng,
                m_state,
                mass,
                thermo_beta
            );
        } else if (init_type == "xyz") {
            //const long effective_frame = override_frame.has_value() ? override_frame.value() : config->init_vel_frame;
            //const long effective_frame = override_frame.value_or(config->init_vel_frame);
            long effective_frame;
            if (override_frame.has_value()) {
                effective_frame = override_frame.value();
            } else {
                if (config->rpmd_config.enabled) {
                    // If RPMD is enabled, set the effective frame to be the first frame of the XYZ file,
                    // and let the RPMD frame selector handle the frame selection for each run later.
                    effective_frame = 0;
                } else {
                    effective_frame = config->init_vel_frame;
                }
            }

            initializer = std::make_unique<XyzMomentumInitializer>(
                filename,
                first_idx,
                unit,
                effective_frame,
                frame_mode,
                m_state,
                mass
            );
        } else {
            throw std::invalid_argument("Unknown momentum initialization method: " + init_type);
        }

        MpiUtils::runCollectively(
            "Momentum initialization",
            [&] {
                initializer->initialize();
            }
        );
    };
    
}

void Simulation::initializePositions(const std::shared_ptr<SimulationConfig>& config)
{
    bindPositionInitFactory(config);
    m_position_init_factory(std::nullopt); // Use the factory with no frame override
}

void Simulation::initializeMomenta(const std::shared_ptr<SimulationConfig>& config)
{
    bindMomentumInitFactory(config);
    m_momentum_init_factory(std::nullopt); // Use the factory with no frame override
}

/*void Simulation::initializePositions(const std::shared_ptr<SimulationConfig>& config)
{
    if (!m_state || !m_rng) {
        throw std::runtime_error("State or RNG is not initialized before position initialization.");
    }

    std::unique_ptr<PositionInitializer> initializer;

    if (config->init_pos_type == "random")
    {
        initializer = std::make_unique<RandomPositionInitializer>(
            m_rng,
            std::shared_ptr<VecArray>(m_state, &m_state->coord),
            m_box_ctx
        );
    }
    else if (config->init_pos_type == "grid")
    {
        initializer = std::make_unique<GridPositionInitializer>(
            std::shared_ptr<VecArray>(m_state, &m_state->coord),
            m_box_ctx
        );
    }
    else if (config->init_pos_type == "xyz")
    {
        initializer = std::make_unique<XyzPositionInitializer>(
            config->init_pos_filename,
            config->this_bead + config->init_pos_index_offset,
            config->init_pos_unit,
            config->init_pos_frame,
            config->init_pos_frame_mode,
            std::shared_ptr<VecArray>(m_state, &m_state->coord),
            m_box_ctx
        );
    }
    else
    {
        throw std::invalid_argument("Unknown position initialization method: " + config->init_pos_type);
    }

    MpiUtils::runCollectively(
        "XYZ coordinate initialization",
        [&] {
            initializer->initialize();
        }
    );


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
    else if (config->init_vel_type == "xyz")
    {
        initializer = std::make_unique<XyzMomentumInitializer>(
            config->init_vel_filename,
            config->this_bead + config->init_vel_index_offset,
            config->init_vel_unit,
            config->init_vel_frame,
            config->init_vel_frame_mode,
            m_state,
            config->mass
        );
    }
    else
    {
        throw std::invalid_argument("Unknown momentum initialization method: " + config->init_vel_type);
    }

    initializer->initialize();
}*/

void Simulation::initializeForces(
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<ForceManager>& force_mgr,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx
)
{
    /*
    if (!m_force_mgr)
    {
        throw std::runtime_error("Force manager must be initialized before force initialization.");
    }

    m_force_mgr->updateForces(*m_state, m_spring_ctx, m_box_ctx);
    */

    force_mgr->updateForces(*state, spring_ctx, box_ctx);
}

void Simulation::setStep(long step)
{
    m_is_thermalization_phase = (step < m_threshold);
}

/*
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
*/

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

void Simulation::initializeObservablesLogger(const std::filesystem::path& filename) {
    // If logger is not yet initialized, create it. Otherwise, reopen the file with the new filename.
    if (!m_obs_logger)
    {
        m_obs_logger = std::make_unique<ObservablesLogger>(
            filename,
            m_bead_ctx.this_bead,
            m_sfreq,
            m_observables
        );
    }
    else
    {
        // Logger already exists, so we just reopen the file with the new filename
        m_obs_logger->reopenFile(filename);
    }
}

void Simulation::initializeDumps(const std::filesystem::path& base_folder)
{
    for (const auto& dump : m_dumps)
    {
        dump->reopenFile(base_folder);
    }
}

void Simulation::calculateObservables() const {
    for (const auto& observable : m_observables) {
        observable->calculate();
    }
}

void Simulation::calculateAndLogObservables(long step) const {
    if (m_obs_logger->isLoggingStep(step)) {
        calculateObservables();
        m_obs_logger->log(step);
        // Resetting the observables is unnecessary since we don't accumulate them across multiple steps.
        // Should one change this, the resetting method can be placed here.
        // resetObservables(); 
    }
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
    if (!m_is_thermalization_phase) {
        dumpStepInfo(step);
        //calculateAndLogObservables(step); <--- TODO: Should be here
    }

    performMolecularDynamicsStep();

    // TODO: Temporarily placed here, to not break the old tests
    if (!m_is_thermalization_phase) {
        calculateAndLogObservables(step);
    }
}

void Simulation::run()
{
    if (m_is_rpmd) {
        runRingPolymerMolecularDynamics();
    } else {
        runPathIntegralMolecularDynamics();
    }
}

void Simulation::runPathIntegralMolecularDynamics() {
    printStatus("Running the simulation", m_bead_ctx.this_bead);
    const double start_time = getWallTime();

    // Initialize the observables logger for the main simulation output
    initializeObservablesLogger(Output::getPimdFilename());
    initializeDumps(Output::FOLDER_NAME);

    // Main loop performing molecular dynamics steps
    for (long step = 0; step <= m_steps; ++step) {
        setStep(step);
        executeStep(step);
    }

    finalizeSimulation(start_time);
}

void Simulation::runRingPolymerMolecularDynamics() {
    if (!m_is_rpmd)
    {
        throw std::runtime_error("Attempted to run RPMD simulation when it is not enabled in config.");
    }

    if (!m_rpmd_frame_selector)
    {
        throw std::runtime_error("RPMD frame selector is not initialized. Ensure that it is created when RPMD is enabled.");
    }

    printStatus("Initializing RPMD", m_bead_ctx.this_bead);
    printStatus(
        std::format("Running {} independent NVE simulations", m_rpmd_context.num_runs),
        m_bead_ctx.this_bead
    );

    const double start_time = getWallTime();

    // Prepare the exact XYZ frame indices based on the number of runs and the fraction of NVT points that need to be discarded
    const std::vector<long>& rpmd_frames = m_rpmd_frame_selector->getRpmdFrameIndices();
    long last_frame_idx = m_rpmd_frame_selector->getNumFrames() - 1;

    for (int run = 0; run < m_rpmd_context.num_runs; ++run) {
        printStatus(std::format("Starting NVE run {}/{} (frame index {}/{})", run + 1, m_rpmd_context.num_runs, rpmd_frames[run], last_frame_idx), m_bead_ctx.this_bead);

        const long frame = rpmd_frames[run];

        // Re-initialize positions and momenta for the current RPMD run
        m_position_init_factory(frame);
        m_momentum_init_factory(frame);

        /// TODO: Re-initialize the forces as well (once we do this in the constructor & update the tests)

        const std::filesystem::path output_filename = Output::getRpmdFilename(run);
        const std::filesystem::path output_folder = Output::getRpmdFolder(run);
        initializeObservablesLogger(output_filename);
        initializeDumps(output_folder);

        // Main loop performing molecular dynamics steps for the current RPMD run
        for (long step = 0; step <= m_steps; ++step) {
            setStep(step);
            executeStep(step);
        }

        printStatus(std::format("Completed NVE run {}/{}", run + 1, m_rpmd_context.num_runs), m_bead_ctx.this_bead);
    }

    const double wall_time = getWallTime() - start_time;

    printStatus(
        std::format(
            "RPMD finished running successfully (Runtime = {:.3} sec)",
            wall_time
        ),
        m_bead_ctx.this_bead
    );

    m_report->writeReport(wall_time); /// TODO: For RPMD we need a separate report
}