#include "simulation.h"
#include "params.h"
#include "core/system_state.h"
#include "core/force_manager.h"
#include "core/random_generators.h"
#include "momentum_initializers.h"
#include "position_initializers.h"
#include "propagators.h"
#include "initializers/observable_initializer.h"
#include "initializers/dump_initializer.h"
#include "thermostats.h"
#include "observables_logger.h"
#include "output_paths.h"

#include <ranges>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <array>
#include <cassert>


#if FACTORIAL_BOSONIC_ALGORITHM
#include "bosonic_exchange/factorial_bosonic_exchange.h"
#else
#include "bosonic_exchange/quadratic_bosonic_exchange.h"
#endif

Simulation::Simulation(int rank, int nproc, const std::string& config_filename): m_step(0)
{
    // Parse the simulation parameters from the configuration (input) file
    const Params params(config_filename, rank);
    m_config = params.load();

    if (!m_config)
    {
        // CR: Doesn't Params throw exceptions? What's the difference between that and returning null?
        throw std::runtime_error("Failed to load configuration");
    }

    // Initialize basic contexts
    // CR: Make sure that they are all disjoint
    const auto thermal_ctx = ThermalContext{
        .beta = m_config->beta,
        .thermo_beta = m_config->thermo_beta
    };

    const auto spring_ctx = SpringContext{
        .omega_p = m_config->omega_p,
        .spring_constant = m_config->spring_constant,
        .beta_half_k = m_config->beta_half_k
    };

    const auto box_ctx = BoxContext{
        .box_size = m_config->box_size,
        .pbc = m_config->pbc
    };

    const auto bead_ctx = BeadContext{
        .nbeads = m_config->nbeads,
        .natoms = m_config->natoms,
        .this_bead = m_config->this_bead
    };

    // CR: What is meaning of the variable prefix m_?
    m_rng = std::make_shared<RandomGenerators>(m_config->seed + rank);
    m_state = std::make_shared<SystemState>(rank, nproc, m_config->natoms, m_config->nbeads);

    initializePositions(m_config, box_ctx, m_state, m_rng);
    initializeMomenta(m_config, m_state, m_rng);

    m_bosonic_exchange = initializeExchange(
        m_config->bosonic,
        thermal_ctx,
        spring_ctx,
        box_ctx,
        bead_ctx,
        m_state
    );

    // Initialize other resources
    m_force_mgr = std::make_shared<ForceManager>(
        m_config,
        m_bosonic_exchange,
        bead_ctx
    );    

    m_normal_modes = initializeNormalModes(m_config, m_state);
    const auto nm_ctx = NormalModesContext{
        .normal_modes = m_normal_modes,
        .couple_to_nm = std::get<bool>(m_config->thermostat_params["nmthermostat"])
    };

    // CR: It's OK imo to pass all of config to the local function (of the class Simulation) that initializes the propagator.
    // CR: Would improve readability of this function.
    // CR: The important bit is that Config & Simulation are not exposed to *other* classes.
    m_propagator = initializePropagator(
        m_config->propagator_type,
        m_config->mass,
        m_config->dt,
        spring_ctx,
        box_ctx,
        bead_ctx,
        m_config->bosonic,
        m_state,
        m_normal_modes,
        m_force_mgr,
        m_bosonic_exchange
    );

    m_thermostat = initializeThermostat(
        thermal_ctx,
        nm_ctx,
        m_config,
        m_state,
        m_rng
    );

    const auto thermostat_ctx = ThermostatContext{
        .thermostat = m_thermostat,
        .thermostat_type = m_config->thermostat_type
    };

    const auto velocity_ctx = VelocityContext{
    .momenta = std::shared_ptr<dVec>(m_state, &m_state->momenta),
    .mass = m_config->mass
    };

    // CR: I'd write this->m_observables for emphasis
    // CR: I think it's OK for these classes to accept config too, because they have
    // CR: no other logic except bridging the gap between config and the simulation.
    // CR: They don't continue to live after initialization
    m_observables = ObservableInitializer(
        m_config,
        m_state,
        m_bosonic_exchange,
        m_force_mgr,
        m_thermostat,
        bead_ctx,
        thermal_ctx,
        spring_ctx,
        box_ctx,
        velocity_ctx,
        thermostat_ctx
    ).createObservables();

    m_dumps = DumpInitializer(m_config, m_state, velocity_ctx).createDumps();
}

Simulation::~Simulation() = default;

std::shared_ptr<BosonicExchangeBase> Simulation::initializeExchange(
    bool bosonic,
    const ThermalContext& thermal_ctx,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    const BeadContext& bead_ctx,
    const std::shared_ptr<SystemState>& state)
{
    bool is_bosonic_bead = bosonic && (bead_ctx.this_bead == 0 || bead_ctx.this_bead == bead_ctx.nbeads - 1);

    // If this is not a bosonic bead, we don't need to initialize the exchange
    if (!is_bosonic_bead) {
        // CR: Can this entire logic be encapsulated in StatisticsManager, which becomes
        // CR: a dynamic rather than static entity?
        // CR: Ideally, the word "bosons" doesn't appear anywhere outside the statistics manager
        // CR: It can be a singleton for ease of access, since its pervasive
        return { nullptr };  /// TODO: Not ideal. Fix this later
    }

    std::shared_ptr<dVec> x_first_bead;
    std::shared_ptr<dVec> x_last_bead;

    if (state->currentBead() == 0)
    {
        // At the first imaginary time slice, the last ("P") slice is the previous one
        x_first_bead = std::shared_ptr<dVec>(state, &state->coord);
        x_last_bead = std::shared_ptr<dVec>(state, &state->prev_coord);
    }
    else
    {
        // At the last imaginary time slice ("P"), the first slice is the next one
        x_first_bead = std::shared_ptr<dVec>(state, &state->next_coord);
        x_last_bead = std::shared_ptr<dVec>(state, &state->coord);
    }

#if FACTORIAL_BOSONIC_ALGORITHM
    return std::make_unique<FactorialBosonicExchange>(
        x_first_bead,
        x_last_bead,
        thermal_ctx,
        spring_ctx,
        box_ctx,
        bead_ctx
    );
#else
    return std::make_shared<BosonicExchange>(
        x_first_bead,
        x_last_bead,
        thermal_ctx,
        spring_ctx,
        box_ctx,
        bead_ctx
    );
#endif
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
    const std::string& propagator_type,
    double mass,
    double dt,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    const BeadContext& bead_ctx,
    bool bosonic,
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<NormalModes>& normal_modes,
    const std::shared_ptr<ForceManager>& force_mgr,
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange)
{
    // JH: propagator_type, mass and dt are passed as arguments instead of the entire config
    //     because I'm confident that they are necessary and sufficient (along with the other arguments)
    //     to initialize any current/future propagator. But perhaps passing the entire config is okay too?
    //     (See the comment in initializeThermostat for a possible reason)
    if (propagator_type == "cartesian")
    {
        return std::make_shared<VelocityVerletPropagator>(
            state,
            force_mgr,
            box_ctx,
            spring_ctx,
            mass,
            dt
        );
    }

    if (propagator_type == "normal_modes")
    {
        return std::make_shared<NormalModesPropagator>(
            state,
            force_mgr,
            box_ctx,
            spring_ctx,
            mass,
            dt,
            bosonic_exchange,
            bead_ctx,
            bosonic,
            normal_modes
        );
    }

    return {nullptr};
}

std::shared_ptr<Thermostat> Simulation::initializeThermostat(
    const ThermalContext& thermal_ctx,
    const NormalModesContext& nm_ctx,
    const std::shared_ptr<SimulationConfig>& config,
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<RandomGenerators>& rng)
{
    // JH: Here I hesitated a little bit to pass individual parameters instead of "config"
    //     because it could be that a new thermostat implementation might need new parameters.
    //     If I had followed the approach of "initializeExchangeState", I would have had to pass
    //     parameters like, e.g., "nchains", which is too specific (it exposes the existence of
    //     Nose-Hoover and its unique parameters). Also, passing all the individual parameters
    //     of all the thermostats would greatly increase the number of arguments and possibly
    //     hinder readability.
    if (config->thermostat_type == "langevin")
    {
        double gamma = std::get<double>(config->thermostat_params["gamma"]);

        return std::make_shared<LangevinThermostat>(
            thermal_ctx,
            nm_ctx,
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
        return std::make_shared<NoseHooverThermostat>(thermal_ctx, nm_ctx, state, nchains, config->dt, config->mass);
    }

    if (config->thermostat_type == "nose_hoover_np")
    {
        return std::make_shared<NoseHooverNpThermostat>(thermal_ctx, nm_ctx, state, nchains, config->dt, config->mass);
    }

    if (config->thermostat_type == "nose_hoover_np_dim")
    {
        return std::make_shared<NoseHooverNpDimThermostat>(thermal_ctx, nm_ctx, state, nchains, config->dt,
                                                           config->mass);
    }

    if (config->thermostat_type == "none")
    {
        return std::make_shared<Thermostat>(thermal_ctx, nm_ctx, state);
    }

    return {nullptr};
}

void Simulation::initializePositions(
    const std::shared_ptr<SimulationConfig>& config,
    // CR: Hmm. I like how the initializer depend on the "contexts" and not the entire config
    // CR: But it's strange to pass here config *and* other things that are included in it
    const BoxContext& box_ctx,
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<RandomGenerators>& rng)
{
    std::unique_ptr<PositionInitializer> initializer;

    if (config->init_pos_type == "random")
    {
        initializer = std::make_unique<RandomPositionInitializer>(
            rng,
            std::shared_ptr<dVec>(state, &state->coord),
            box_ctx
        );
    }
    else if (config->init_pos_type == "grid")
    {
        initializer = std::make_unique<GridPositionInitializer>(
            std::shared_ptr<dVec>(state, &state->coord),
            box_ctx
        );
    }
    else if (config->init_pos_type == "xyz")
    {
        initializer = std::make_unique<XyzPositionInitializer>(
            config->init_pos_filename,
            config->this_bead + config->init_pos_index_offset,
            std::shared_ptr<dVec>(state, &state->coord),
            box_ctx
        );
    }
    else
    {
        throw std::invalid_argument("Unknown position initialization method: " + config->init_pos_type);
    }

    initializer->initialize();

    // Communicate the new coordinates to the neighboring processes
    // CR: Good comment
    state->updateNeighboringCoordinates();
}

void Simulation::initializeMomenta(
    const std::shared_ptr<SimulationConfig>& config,
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<RandomGenerators>& rng)
{
    std::unique_ptr<MomentumInitializer> initializer;

    if (config->init_vel_type == "random")
    {
        initializer = std::make_unique<MaxwellBoltzmannMomentumInitializer>(
            rng,
            state,
            config->mass,
            config->thermo_beta
        );
    }
    else if (config->init_vel_type == "manual")
    {
        initializer = std::make_unique<ManualMomentumInitializer>(
            config->init_vel_filename,
            config->init_vel_index_offset,
            state,
            config->mass
        );
    }
    else
    {
        throw std::invalid_argument("Unknown momentum initialization method: " + config->init_vel_type);
    }

    initializer->initialize();
}

int Simulation::getStep() const
{
    return m_step;
}

void Simulation::setStep(const int step)
{
    m_step = step;
}

void Simulation::run()
{
    printStatus("Running the simulation", m_config->this_bead);

    MPI_Barrier(MPI_COMM_WORLD);
    const double sim_exec_time_start = MPI_Wtime();

    std::filesystem::create_directory(Output::FOLDER_NAME);
    // Initialize the output file for the observables
    ObservablesLogger obs_logger(Output::MAIN_FILENAME, m_config->this_bead, m_observables);

    // Initialize the files for the dumps (e.g., xyz, dat)
    // CR: extract method, initializeDumps()
    // CR: Also, why aren't they initialized with everything else
    for (const auto& dump : m_dumps)
    {
        dump->initialize();
    }

    // Main loop performing molecular dynamics steps
    for (long step = 0; step <= m_config->steps; ++step)
    {
        // CR: extract method - step()

        // CR: Is it possible to separate the things that change each step to a different object?
        setStep(step);

        // Reset the observables at the beginning of each step
        /// TODO: Do we need this? What if we want accumulation? Does it account for frequency?
        for (const auto& observable : m_observables)
        {
            observable->resetValues();
        }

        // Dump the desired quantities (e.g., coordinates, forces, etc.) at the specified frequency
        // CR: Is this comment helpful?
        // CR: extract method dumpStepInfo() or the like
        for (const auto& dump : m_dumps)
        {
            dump->output(step);
        }

        // Perform a thermostat step
        // CR: Is this comment helpful?
        m_thermostat->step();

        // If fixcom=true, the center of mass of the ring polymers is fixed during the simulation
        // CR: Is this comment helpful? etc.
        if (m_config->fixcom)
        {
            m_state->zeroMomentum();
        }

        // CR: is this comment useful?
        // Perform a time propagation step
        m_propagator->step();

        // CR: is this comment useful?
        // Perform a thermostat step
        m_thermostat->step();

        // Zero momentum after every thermostat step (if needed)
        if (m_config->fixcom)
        {
            m_state->zeroMomentum();
        }

        // If we have not reached the thermalization threshold, skip to the next step (thermalization stage)
        if (step < m_config->threshold) {
            // CR: If your function has a continue statement, it must be very short, and the flow very clear
            // CR: Can be solved by moving everything after this to a separate method
            // CR: or: calculateObservableIfAskedAtThisStep
            continue;
        }

        // Calculate the observables (production stage)
        for (const auto& observable : m_observables)
        {
            observable->calculate();
        }

        // Save the observables at the specified frequency
        if (step % m_config->sfreq == 0)
        {
            obs_logger.log(step);
        }
    }

    // CR: What's the purpose of the barrier?
    MPI_Barrier(MPI_COMM_WORLD);

    // CR: extact method
    const double sim_exec_time_end = MPI_Wtime();
    const double wall_time = sim_exec_time_end - sim_exec_time_start;

    printStatus(std::format("Simulation finished running successfully (Runtime = {:.3} sec)", wall_time),
                m_config->this_bead);

    printReport(wall_time);
}

void Simulation::printReport(double wall_time) const
{
    if (m_config->this_bead != 0)
        return;

    std::ofstream report_file;
    report_file.open(std::format("{}/report.txt", Output::FOLDER_NAME), std::ios::out | std::ios::app);

    report_file << "---------\nParameters\n---------\n";

    if (m_config->bosonic)
    {
        report_file << formattedReportLine("Statistics", "Bosonic");
        std::string bosonic_alg_name = "Feldman-Hirshberg";

#if FACTORIAL_BOSONIC_ALGORITHM
        bosonic_alg_name = "Naive";
#endif

        report_file << formattedReportLine("Bosonic algorithm", bosonic_alg_name);
    }
    else
    {
        report_file << formattedReportLine("Statistics", "Boltzmannonic");
    }

    report_file << formattedReportLine("Time propagation algorithm", m_config->propagator_type);
    report_file << formattedReportLine("Periodic boundary conditions", m_config->pbc);
    report_file << formattedReportLine("Dimension", NDIM);
    report_file << formattedReportLine("Seed", m_config->seed);
    report_file << formattedReportLine("Coordinate initialization method", m_config->init_pos_type);
    report_file << formattedReportLine("Momentum initialization method", m_config->init_vel_type);
    report_file << formattedReportLine("Number of atoms", m_config->natoms);
    report_file << formattedReportLine("Number of beads", m_config->nbeads);

    /// TODO: Limit the number of digits in the output
    double out_temperature = Units::convertToUser("temperature", "kelvin", m_config->temperature);
    report_file << formattedReportLine("Temperature", std::format("{} kelvin", out_temperature));

    double out_sys_size = Units::convertToUser("length", "angstrom", m_config->box_size);
    report_file << formattedReportLine("Linear size of the system", std::format("{} angstroms", out_sys_size));

    double out_mass = Units::convertToUser("mass", "dalton", m_config->mass);
    report_file << formattedReportLine("Mass", std::format("{} amu", out_mass));

    report_file << formattedReportLine("Total number of MD steps", m_config->steps);
    report_file << formattedReportLine("Interaction potential name", m_config->int_pot_name);
    report_file << formattedReportLine("External potential name", m_config->ext_pot_name);

    report_file << "---------\nFeatures\n---------\n";
    report_file << formattedReportLine("Minimum image convention", MINIM);
    report_file << formattedReportLine("Wrapping of coordinates", WRAP);
    report_file << formattedReportLine("Using i-Pi convention", IPI_CONVENTION);

    report_file << "---------\n";
    report_file << formattedReportLine("Wall time", std::format("{:%T}",
                                                                std::chrono::duration<double>(wall_time)
                                       ));
    report_file << formattedReportLine("Wall time per step (sec)",
                                       std::format("{:.5e}", wall_time / m_config->steps));

    report_file.close();
}
