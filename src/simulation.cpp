#include "simulation.h"

#include "common.h"
#include "params.h"

#include "core/system_state.h"
#include "core/force_manager.h"
#include "core/random_generators.h"
#include "momentum_initializers.h"
#include "position_initializers.h"

//#include "contexts/spring_context.h"
//#include "contexts/velocity_context.h"
//#include "contexts/thermal_context.h"
//#include "contexts/box_context.h"
//#include "contexts/normal_modes_context.h"

#include "propagators.h"
#include "dumps.h"
#include "observables.h"
#include "thermostats.h"

#include <ranges>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <array>
#include <cassert>

#include "observables_logger.h"

Simulation::Simulation(int rank, int nproc, const std::string& config_filename): m_step(0)
{
    // Load the simulation parameters from the configuration (input) file
    const Params params(config_filename, rank);
    const std::shared_ptr<SimulationConfig> config = params.load();

    if (!config) {
        throw std::runtime_error("Failed to load configuration");
    }

    // Initialize the random number generator (each process has a unique seed)
    const auto rng = std::make_shared<RandomGenerators>(config->seed + rank);

    // Initialize the simulation state
    const auto state = std::make_shared<SystemState>(rank, nproc, config->natoms, config->nbeads);

    // Initialize particle positions and velocities
    initializePositions(config, state, rng);
    initializeMomenta(config, state, rng);

    // Communicate the new coordinates to the neighboring processes
    state->updateNeighboringCoordinates();

    // Initialize the force manager
    const auto force_mgr = std::make_shared<ForceManager>(config);

    // Initialize the exchange state
    const auto exchange_state = initializeExchangeState(config, state);

    // Initialize normal modes (if necessary)
    const auto normal_modes = initializeNormalModes(config, state);

    // Initialize the time propagation scheme
    const auto propagator = initializePropagator(config, state, normal_modes, force_mgr, exchange_state);

    // Initialize the thermostat
    const auto thermostat = initializeThermostat(config, state, normal_modes, rng);

    // Initialize the observables
    const auto observables = initializeObservables(config, state, exchange_state, force_mgr, thermostat);

    // Initialize the dumps
    const auto dumps = initializeDumps(config, state);

    // Load the simulation context with the initialized objects
    m_context = SimulationResources{
        .config = config,
        .state = state,
        .exchange_state = exchange_state,
        .rng = rng,
        .force_mgr = force_mgr,
        .normal_modes = normal_modes,
        .propagator = propagator,
        .thermostat = thermostat,
        .observables = observables,
        .dumps = dumps
    };
}

Simulation::~Simulation() = default;

std::shared_ptr<ExchangeState> Simulation::initializeExchangeState(
    const std::shared_ptr<SimulationConfig>& config,
    const std::shared_ptr<SystemState>& state)
{
    std::shared_ptr<dVec> x_first_bead;
    std::shared_ptr<dVec> x_last_bead;
    std::shared_ptr<dVec> x_neighbor_bead;

    if (config->this_bead == 0)
    {
        x_first_bead = std::shared_ptr<dVec>(state, &state->coord);
        x_last_bead = std::shared_ptr<dVec>(state, &state->prev_coord);
    } else
    {
        x_first_bead = std::shared_ptr<dVec>(state, &state->next_coord);
        x_last_bead = std::shared_ptr<dVec>(state, &state->coord);
    }

    return std::make_shared<ExchangeState>(
        x_first_bead,
        x_last_bead,
        ThermalContext{
            .beta = config->beta,
            .thermo_beta = config->thermo_beta
        },
        SpringContext{
            .omega_p = config->omega_p,
            .spring_constant = config->spring_constant,
            .beta_half_k = config->beta_half_k
        },
        BoxContext{
            .box_size = config->box_size,
            .pbc = config->pbc
        },
        BeadContext{
            .nbeads = config->nbeads,
            .natoms = config->natoms,
            .this_bead = config->this_bead,
        },
        config->bosonic
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
    const std::shared_ptr<ForceManager>& force_mgr,
    const std::shared_ptr<ExchangeState>& exchange_state)
{
    /// TODO: Pass to initializePropagator the spring_context, mass and dt instead of the entire config?
    SpringContext spring_context{
        .omega_p = config->omega_p,
        .spring_constant = config->spring_constant,
        .beta_half_k = config->beta_half_k
    };

    if (config->propagator_type == "cartesian")
    {
        return std::make_shared<VelocityVerletPropagator>(
            state,
            force_mgr,
            exchange_state,
            spring_context,
            config->mass,
            config->dt
        );
    }

    if (config->propagator_type == "normal_modes")
    {
        return std::make_shared<NormalModesPropagator>(
            state,
            force_mgr,
            exchange_state,
            spring_context,
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
    const std::shared_ptr<NormalModes>& normal_modes,
    const std::shared_ptr<RandomGenerators>& rng)
{
    ThermalContext thermal_ctx{
        .beta = config->beta,
        .thermo_beta = config->thermo_beta
    };

    NormalModesContext nm_ctx{
        .normal_modes = normal_modes,
        .couple_to_nm = std::get<bool>(config->thermostat_params["nmthermostat"])
    };

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
        return std::make_shared<NoseHooverNpDimThermostat>(thermal_ctx, nm_ctx, state, nchains, config->dt, config->mass);
    }

    if (config->thermostat_type == "none")
    {
        return std::make_shared<Thermostat>(thermal_ctx, nm_ctx, state);
    }

    return {nullptr};
}

std::vector<std::shared_ptr<Observable>> Simulation::initializeObservables(
    const std::shared_ptr<SimulationConfig>& config,
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<ExchangeState>& exchange_state,
    const std::shared_ptr<ForceManager>& force_mgr,
    const std::shared_ptr<Thermostat>& thermostat)
{
    std::vector<std::shared_ptr<Observable>> observables;
    observables.reserve(config->observables_list.size()); // Pre-allocate for efficiency

    for (const auto& [obs_name, obs_unit] : config->observables_list)
    {
        // Skip disabled observables
        if (obs_unit == "off")
        {
            continue;
        }

        // Determine the correct unit to use
        std::string correct_unit = (obs_unit != "none") ? obs_unit : "";

        // Use enum-based switch for better performance and readability
        enum ObservableType : std::int8_t
        {
            ENERGY,
            CLASSICAL,
            BOSONIC,
            GSF,
            UNKNOWN
        };

        // CR: What is the role of the separation between observable types?
        // JH: They represent different types of observable categories. For example, if
        //      one requests "energy", it makes sense to calculate the different types of energies (primitive KE, virial KE, PE, etc.),
        //      but not calculate other estimators (unless the user asks for it).
        // CR: "classical" seems to need access to most of the things?
        // JH: Everything passed to it is eventually used.

        ObservableType type = UNKNOWN;
        if (obs_name == "energy") type = ENERGY;
        else if (obs_name == "classical") type = CLASSICAL;
        else if (obs_name == "bosonic") type = BOSONIC;
        else if (obs_name == "gsf") type = GSF;

        switch (type)
        {
        case ENERGY:
            observables.push_back(std::make_shared<EnergyObservable>(
                exchange_state,
                std::shared_ptr<const dVec>(state, &state->coord),
                std::shared_ptr<const dVec>(state, &state->prev_coord),
                force_mgr,
                BeadContext{
                    .nbeads = config->nbeads,
                    .natoms = config->natoms,
                    .this_bead = config->this_bead
                },
                ThermalContext{
                    .beta = config->beta,
                    .thermo_beta = config->thermo_beta
                },
                SpringContext{
                    .omega_p = config->omega_p,
                    .spring_constant = config->spring_constant,
                    .beta_half_k = config->beta_half_k
                },
                BoxContext{
                    .box_size = config->box_size,
                    .pbc = config->pbc
                },
                config->sfreq, /// TODO: Generalize the save frequency option to dumps, observables, etc.
                correct_unit
            ));
            break;

        case CLASSICAL:
            observables.push_back(std::make_shared<ClassicalObservable>(
                std::shared_ptr<const dVec>(state, &state->coord),
                std::shared_ptr<const dVec>(state, &state->prev_coord),
                exchange_state,
                VelocityContext{
                    .momenta = std::shared_ptr<const dVec>(state, &state->momenta),
                    .mass = config->mass
                },
                ThermostatContext{
                    .thermostat = thermostat,
                    .thermostat_type = config->thermostat_type
                },
                BeadContext{
                    .nbeads = config->nbeads,
                    .natoms = config->natoms,
                    .this_bead = config->this_bead
                },
                SpringContext{
                    .omega_p = config->omega_p,
                    .spring_constant = config->spring_constant,
                    .beta_half_k = config->beta_half_k
                },
                BoxContext{
                    .box_size = config->box_size,
                    .pbc = config->pbc
                },
                config->sfreq, /// TODO: Generalize the save frequency option to dumps, observables, etc.
                correct_unit
            ));
            break;

        case BOSONIC:
            if (!config->bosonic)
            {
                throw std::runtime_error("Bosonic observables require bosonic simulation mode");
            }
            observables.push_back(std::make_shared<BosonicObservable>(
                exchange_state,
                config->sfreq, /// TODO: Generalize the save frequency option to dumps, observables, etc.
                correct_unit
            ));
            break;

        case GSF:
            observables.push_back(std::make_shared<GSFActionObservable>(
                std::shared_ptr<const dVec>(state, &state->coord),
                force_mgr,
                BeadContext{
                    .nbeads = config->nbeads,
                    .natoms = config->natoms,
                    .this_bead = config->this_bead
                },
                ThermalContext{
                    .beta = config->beta,
                    .thermo_beta = config->thermo_beta
                },
                SpringContext{
                    .omega_p = config->omega_p,
                    .spring_constant = config->spring_constant,
                    .beta_half_k = config->beta_half_k
                },
                config->sfreq, /// TODO: Generalize the save frequency option to dumps, observables, etc.
                correct_unit
            ));
            break;

        case UNKNOWN:
        default:
            throw std::runtime_error("Unknown observable type: " + obs_name);
        }
    }

    return observables;
}

std::vector<std::shared_ptr<Dump>> Simulation::initializeDumps(
    const std::shared_ptr<SimulationConfig>& config,
    const std::shared_ptr<SystemState>& state)
{
    std::vector<std::shared_ptr<Dump>> dumps;
    dumps.reserve(config->dumps_list.size()); // Pre-allocate for efficiency

    // Initialize the dump files based on the configuration
    for (const auto& [dump_name, dump_unit] : config->dumps_list)
    {
        if (dump_unit == "off" || dump_unit == "false")
        {
            continue;
        }

        // Check if using default unit
        const bool use_default_unit = (dump_unit == "none" || dump_unit == "true" || dump_unit == "on");

        // Determine the correct unit to use
        //std::string correct_unit = (dump_unit != "none" && dump_unit != "true" && dump_unit != "on") ? dump_unit : "atomic_unit";
        std::string correct_unit = use_default_unit ? "atomic_unit" : dump_unit;

        // Use enum-based switch for better performance and readability
        enum DumpType : std::int8_t
        {
            POSITION,
            VELOCITY,
            FORCE,
            UNKNOWN
        };

        DumpType type = UNKNOWN;
        if (dump_name == "positions") type = POSITION;
        else if (dump_name == "velocities") type = VELOCITY;
        else if (dump_name == "forces") type = FORCE;

        switch (type)
        {
        case POSITION:
            dumps.push_back(std::make_shared<PositionDump>(
                std::shared_ptr<dVec>(state, &state->coord),
                config->this_bead,
                config->sfreq, /// TODO: Generalize the save frequency option to dumps, observables, etc.
                correct_unit
            ));
            break;

        case VELOCITY:
            dumps.push_back(std::make_shared<VelocityDump>(
                // CR: this struct seems to be used only in velocity_dump.cpp?
                // JH: Yes, it is used only in velocity_dump.cpp. It is a context for the VelocityDump class.
                VelocityContext{
                    .momenta = std::shared_ptr<dVec>(state, &state->momenta),
                    .mass = config->mass
                },
                config->this_bead,
                config->sfreq, /// TODO: Generalize the save frequency option to dumps, observables, etc.
                correct_unit
            ));
            break;

        case FORCE:
            dumps.push_back(std::make_shared<ForceDump>(
                state,
                config->this_bead,
                config->sfreq, /// TODO: Generalize the save frequency option to dumps, observables, etc.
                correct_unit
            ));
            break;

        case UNKNOWN:
        default:
            throw std::runtime_error("Unknown dump type: " + dump_name);
        }
    }

    return dumps;
}

void Simulation::initializePositions(
    const std::shared_ptr<SimulationConfig>& config,
    const std::shared_ptr<SystemState>& state,
    const std::shared_ptr<RandomGenerators>& rng)
{
    std::unique_ptr<PositionInitializer> initializer;

    if (config->init_pos_type == "random")
    {
        initializer = std::make_unique<RandomPositionInitializer>(
            rng,
            std::shared_ptr<dVec>(state, &state->coord),
            config->box_size
        );
    }
    else if (config->init_pos_type == "grid")
    {
        initializer = std::make_unique<GridPositionInitializer>(
            std::shared_ptr<dVec>(state, &state->coord),
            config->box_size
        );
    }
    else if (config->init_pos_type == "xyz")
    {
        initializer = std::make_unique<XyzPositionInitializer>(
            config->init_pos_filename,
            config->this_bead + config->init_pos_index_offset,
            std::shared_ptr<dVec>(state, &state->coord),
            config->box_size
        );
    }
    else
    {
        throw std::invalid_argument("Unknown position initialization method: " + config->init_pos_type);
    }

    initializer->initialize();

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
    printStatus("Running the simulation", m_context.config->this_bead);

    MPI_Barrier(MPI_COMM_WORLD);
    const double sim_exec_time_start = MPI_Wtime();

    std::filesystem::create_directory(Output::FOLDER_NAME);
    // Initialize the output file for the observables
    ObservablesLogger obs_logger(Output::MAIN_FILENAME, m_context.config->this_bead, m_context.observables);

    // Initialize the files for the dumps (e.g., xyz, dat)
    for (const auto& dump : m_context.dumps)
    {
        dump->initialize();
    }

    // Main loop performing molecular dynamics steps
    for (long step = 0; step <= m_context.config->steps; ++step)
    {
        setStep(step);

        // Reset the observables at the beginning of each step
        /// TODO: Do we need this? What if we want accumulation? Does it account for frequency?
        for (const auto& observable : m_context.observables)
        {
            observable->resetValues();
        }

        // Dump the desired quantities (e.g., coordinates, forces, etc.) at the specified frequency
        for (const auto& dump : m_context.dumps)
        {
            dump->output(step);
        }

        // Perform a thermostat step
        m_context.thermostat->step();

        // If fixcom=true, the center of mass of the ring polymers is fixed during the simulation
        if (m_context.config->fixcom)
        {
            m_context.state->zeroMomentum();
        }

        // Perform a time propagation step
        m_context.propagator->step();

        // Perform a thermostat step
        m_context.thermostat->step();

        // Zero momentum after every thermostat step (if needed)
        if (m_context.config->fixcom)
        {
            m_context.state->zeroMomentum();
        }

#if PROGRESS
        printProgress(step, steps, m_context.config->this_bead);
#endif

        // If we have not reached the thermalization threshold, skip to the next step (thermalization stage)
        if (step < m_context.config->threshold)
        {
            continue;
        }

        // Calculate the observables (production stage)
        for (const auto& observable : m_context.observables)
        {
            observable->calculate();
        }

        // Save the observables at the specified frequency
        if (step % m_context.config->sfreq == 0)
        {
            obs_logger.log(step);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    const double sim_exec_time_end = MPI_Wtime();
    const double wall_time = sim_exec_time_end - sim_exec_time_start;

    printStatus(std::format("Simulation finished running successfully (Runtime = {:.3} sec)", wall_time),
                m_context.config->this_bead);

    printReport(wall_time);
}

void Simulation::printReport(double wall_time) const
{
    if (m_context.config->this_bead != 0)
        return;

    std::ofstream report_file;
    report_file.open(std::format("{}/report.txt", Output::FOLDER_NAME), std::ios::out | std::ios::app);

    report_file << "---------\nParameters\n---------\n";

    if (m_context.config->bosonic)
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

    report_file << formattedReportLine("Time propagation algorithm", m_context.config->propagator_type);
    report_file << formattedReportLine("Periodic boundary conditions", m_context.config->pbc);
    report_file << formattedReportLine("Dimension", NDIM);
    report_file << formattedReportLine("Seed", m_context.config->seed);
    report_file << formattedReportLine("Coordinate initialization method", m_context.config->init_pos_type);
    report_file << formattedReportLine("Momentum initialization method", m_context.config->init_vel_type);
    report_file << formattedReportLine("Number of atoms", m_context.config->natoms);
    report_file << formattedReportLine("Number of beads", m_context.config->nbeads);

    /// TODO: Limit the number of digits in the output
    double out_temperature = Units::convertToUser("temperature", "kelvin", m_context.config->temperature);
    report_file << formattedReportLine("Temperature", std::format("{} kelvin", out_temperature));

    double out_sys_size = Units::convertToUser("length", "angstrom", m_context.config->box_size);
    report_file << formattedReportLine("Linear size of the system", std::format("{} angstroms", out_sys_size));

    double out_mass = Units::convertToUser("mass", "dalton", m_context.config->mass);
    report_file << formattedReportLine("Mass", std::format("{} amu", out_mass));

    report_file << formattedReportLine("Total number of MD steps", m_context.config->steps);
    report_file << formattedReportLine("Interaction potential name", m_context.config->int_pot_name);
    report_file << formattedReportLine("External potential name", m_context.config->ext_pot_name);

    report_file << "---------\nFeatures\n---------\n";
    report_file << formattedReportLine("Minimum image convention", MINIM);
    report_file << formattedReportLine("Wrapping of coordinates", WRAP);
    report_file << formattedReportLine("Using i-Pi convention", IPI_CONVENTION);

    report_file << "---------\n";
    report_file << formattedReportLine("Wall time", std::format("{:%T}",
                                                                std::chrono::duration<double>(wall_time)
                                       ));
    report_file << formattedReportLine("Wall time per step (sec)",
                                       std::format("{:.5e}", wall_time / m_context.config->steps));

    report_file.close();
}