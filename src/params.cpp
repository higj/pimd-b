#include "params.h"
#include "../libs/string_utils.h"

#include <regex>
#include <format>
#include <filesystem>

Params::Params(const std::string& filename, const int& rank) : m_reader(filename), m_rank(rank) {
    printStatus("Initializing the simulation parameters", rank);

    if (m_reader.ParseError() < 0)
        throw std::invalid_argument(std::format("Unable to read the configuration file {}", filename));
}

std::shared_ptr<SimulationConfig> Params::load() const {
    auto config = std::make_shared<SimulationConfig>();

    config->this_bead = m_rank;

    loadSimulationParams(*config);
    loadSystemParams(*config);
    loadPropagatorParams(*config);
    loadThermostatParams(*config);
    loadCoordInitParams(*config);
    loadMomentaInitParams(*config);
    loadExternalPotentialParams(*config);
    loadInteractionPotentialParams(*config);
    loadOutputParams(*config);
    loadObservableParams(*config);

    return config;
}

/**
 * Load basic parameters pertaining to the simulation.
 *
 * @param config Config object to load parameters into.
 */
void Params::loadSimulationParams(SimulationConfig& config) const {
    config.dt = Units::getQuantity("time", m_reader.Get(Sections::SIMULATION, "dt", "1.0 femtosecond"));
    if (config.dt <= 0.0)
        throw std::invalid_argument(std::format("Invalid time-step ({} is not positive)", config.dt));

    config.threshold = m_reader.GetReal(Sections::SIMULATION, "threshold", 0.0);
    if (config.threshold < 0)
        throw std::invalid_argument(std::format("Invalid threshold ({} is negative)", config.threshold));

    config.steps = static_cast<long>(std::stod(m_reader.Get(Sections::SIMULATION, "steps", "1e5")));
    if (config.steps < 1)
        throw std::invalid_argument(std::format("Invalid number of steps ({}<1)", config.steps));

    config.threshold *= config.steps;  // Threshold in terms of steps

    //config.sfreq = m_reader.GetLong(Sections::SIMULATION, "sfreq", 1000); /// @todo: Add support for scientific notation
    config.sfreq = static_cast<long>(std::stod(m_reader.Get(Sections::SIMULATION, "sfreq", "1e3")));
    if (config.sfreq < 1)
        throw std::invalid_argument(std::format("Invalid value for observable recording frequency ({}<1)", config.sfreq));

    config.nbeads = m_reader.GetInteger(Sections::SIMULATION, "nbeads", 4);
    if (config.nbeads < 1)
        throw std::invalid_argument(std::format("Invalid number of beads ({}<1)", config.nbeads));

    // unsigned int seed = static_cast<unsigned int>(time(nullptr))
    config.seed = static_cast<unsigned int>(std::stod(m_reader.Get(Sections::SIMULATION, "seed", "1234")));

    // Bosonic or distinguishable simulation?
    config.bosonic = m_reader.GetBoolean(Sections::SIMULATION, "bosonic", false);
    // Fix the center of mass?
    config.fixcom = m_reader.GetBoolean(Sections::SIMULATION, "fixcom", true);
    // Enable periodic boundary conditions?
    config.pbc = m_reader.GetBoolean(Sections::SIMULATION, "pbc", false);
}

/**
 * Load basic parameters pertaining to the classical ring-polymer system.
 *
 * @param config Config object to load parameters into.
 */
void Params::loadSystemParams(SimulationConfig& config) const {
    config.natoms = m_reader.GetInteger(Sections::SYSTEM, "natoms", 1);
    if (config.natoms < 1)
        throw std::invalid_argument(std::format("Invalid number of particles (specified {})", config.natoms));

    config.mass = Units::getQuantity("mass", m_reader.Get(Sections::SYSTEM, "mass", "1.0 dalton"));
    if (config.mass <= 0.0)
        throw std::invalid_argument(std::format("The provided mass ({0:4.3f}) is unphysical!", config.mass));

    config.box_size = Units::getQuantity("length", m_reader.Get(Sections::SYSTEM, "size", "1.0 picometer"));
    if (config.box_size <= 0.0)
        throw std::invalid_argument(std::format("The provided system size ({0:4.3f}) is unphysical!", config.box_size));
}

/**
 * Load basic parameters pertaining to the time propagation scheme.
 *
 * @param config Config object to load parameters into.
 */
void Params::loadPropagatorParams(SimulationConfig& config) const {
    // Implemented time propagators:
    //  "cartesian": regular velocity Verlet algorithm, propagating the plain Cartesian coordinates
    //  "normal_modes": a velocity Verlet algorithm that propagates the normal modes
    const StringsList allowed_propagators = { "cartesian", "normal_modes" };
    config.propagator_type = m_reader.GetString(Sections::SIMULATION, "propagator", "cartesian");
    // TODO: Slight danger that config.bosonic might not be properly defined at this point
    if (config.bosonic && config.propagator_type == "normal_modes")
        throw std::invalid_argument("Normal modes propagation is currently not available for bosons!");

    if (!StringUtils::labelInArray(config.propagator_type, allowed_propagators))
        throw std::invalid_argument(std::format("The specified time propagator ({}) is not supported!", config.propagator_type));
}

/**
 * Load parameters pertaining to the thermostat.
 *
 * @param config Config object to load parameters into.
 */
void Params::loadThermostatParams(SimulationConfig& config) const {
    // Setup allowed thermostats with documentation
    const StringsList allowed_thermostats = {
        "langevin",            // A Langevin thermostat coupled to the Cartesian coordinates
        "nose_hoover",         // A single Nose-Hoover chain coupled to the whole system
        "nose_hoover_np",      // A unique Nose-Hoover chain coupled to each particle
        "nose_hoover_np_dim",  // A unique Nose-Hoover chain coupled to each Cartesian coordinate of each particle
        "none"                 // No thermostat (NVE simulation)
    };

    // Read and validate thermostat type first (this drives other validations)
    config.thermostat_type = m_reader.GetString(Sections::SIMULATION, "thermostat", "error");
    if (config.thermostat_type == "error") {
        throw std::invalid_argument("Thermostat must be specified!");
    }
    if (!StringUtils::labelInArray(config.thermostat_type, allowed_thermostats)) {
        throw std::invalid_argument(std::format("Unsupported thermostat type ({})", config.thermostat_type));
    }

    // Validate temperature. When no thermostat is used, temperature value may still be used for initializing velocities
    config.temperature = Units::getQuantity("temperature", m_reader.Get(Sections::SYSTEM, "temperature", "1.0 kelvin"));
    if (config.temperature <= 0) {
        throw std::invalid_argument(std::format("The specified temperature ({0:4.3f} kelvin) is not positive", config.temperature));
    }

    // Define additional parameters based on rudimentary config parameters
    config.beta = 1.0 / (Constants::kB * config.temperature);

#if IPI_CONVENTION
    // i-Pi convention [J. Chem. Phys. 133, 124104 (2010)]
    config.thermo_beta = config.beta / config.nbeads;
    config.omega_p = config.nbeads / (config.beta * Constants::hbar);
#else
    // Tuckerman convention
    config.thermo_beta = config.beta;
    config.omega_p = sqrt(config.nbeads) / (config.beta * Constants::hbar);
#endif

    config.spring_constant = config.mass * config.omega_p * config.omega_p;
    config.beta_half_k = config.thermo_beta * 0.5 * config.spring_constant;

    // Determine if this is a Nose-Hoover type thermostat
    const bool is_nose_hoover = config.thermostat_type.find("nose_hoover") != std::string::npos;

    // Handle nchains parameter (only for Nose-Hoover thermostats)
    if (m_reader.HasValue(Sections::SIMULATION, "nchains") && !is_nose_hoover) {
        throw std::invalid_argument("nchains can only be used with Nose-Hoover thermostats!");
    }

    // Only set nchains for Nose-Hoover thermostats
    if (is_nose_hoover) {
        int nchains = m_reader.GetInteger(Sections::SIMULATION, "nchains", 4);
        if (nchains < 1) {
            throw std::invalid_argument(std::format("Invalid number of Nose-Hoover chains ({}<1)", nchains));
        }
        config.thermostat_params["nchains"] = nchains;
    }

    // Handle normal mode thermostat coupling
    bool nmthermostat = m_reader.GetBoolean(Sections::SIMULATION, "nmthermostat", false);
    if (nmthermostat && config.thermostat_type == "none") {
        throw std::invalid_argument("nmthermostat cannot be used in nve ensemble!");
    }

    // Only set nmthermostat if we're using a thermostat
    if (config.thermostat_type != "none") {
        config.thermostat_params["nmthermostat"] = nmthermostat;
    }

    // Handle gamma parameter (only relevant for Langevin thermostat)
    if (config.thermostat_type == "langevin") {
        const bool gamma_specified = m_reader.HasValue(Sections::SIMULATION, "gamma");
        double gamma;

        if (gamma_specified) {
            // User specified gamma, validate it's positive
            gamma = m_reader.GetReal(Sections::SIMULATION, "gamma", 0.0);
            if (gamma <= 0.0) {
                throw std::invalid_argument(std::format("Invalid gamma value ({0:g}), must be positive", gamma));
            }
        } else {
            // Set default for Langevin thermostat when not specified by user
            gamma = 1 / (100.0 * config.dt);
        }

        config.thermostat_params["gamma"] = gamma;
    }
}

/**
 * Load parameters pertaining to the coordinate initialization.
 *
 * @param config Config object to load parameters into.
 */
void Params::loadCoordInitParams(SimulationConfig& config) const {
    std::string init_pos_type, init_pos_specification;

    if (!StringUtils::parseTokenParentheses(
        m_reader.Get(Sections::SIMULATION, "initial_position", "random"), 
        init_pos_type,
        init_pos_specification
        )) {
        throw std::invalid_argument("Invalid coordinate initialization method");
    }

    const StringsList allowed_coord_init_methods = { "random", "xyz", "grid" }; /// TODO: Use cell instead of grid?

    if (!StringUtils::labelInArray(init_pos_type, allowed_coord_init_methods))
        throw std::invalid_argument(std::format("The specified coordinate initialization method ({}) is not supported!",
            init_pos_type));

    if (init_pos_type == "xyz") {
        try {
            // If the user provided a filename format, try to format it using the first bead index
            const int dummy = 0;
            std::string formatted_filename = std::vformat(init_pos_specification, std::make_format_args(dummy));
            /// TODO: What if both 0 and 1 don't exist?
            config.init_pos_index_offset = static_cast<int>(!std::filesystem::exists(formatted_filename));
        } catch (const std::format_error&) {
            throw std::invalid_argument(
                std::format("The filename format ({}) for coordinate initialization is invalid!",
                    init_pos_specification)
            );
        } catch (...) {
            throw std::runtime_error(
                std::format("Filename format ({}) for coordinate initialization validation failed",
                    init_pos_specification)
            );
        }
    }

    config.init_pos_type = init_pos_type;

    if (init_pos_type == "xyz")
        config.init_pos_filename = init_pos_specification;
}

/**
 * Load parameters pertaining to the momentum initialization.
 *
 * @param config Config object to load parameters into.
 */
void Params::loadMomentaInitParams(SimulationConfig& config) const {
    // Allowed velocity initialization methods:
    //  "random": samples from Maxwell-Boltzmann distribution
    //  "manual": reads from vel_X.dat files
    //  "manual(format)": reads from format(X) files
    std::string init_vel_type, init_vel_specification;

    if (!StringUtils::parseTokenParentheses(
        m_reader.Get(Sections::SIMULATION, "initial_velocity", "random"), 
        init_vel_type,
        init_vel_specification
        )) {
        throw std::invalid_argument("Invalid velocity initialization method");
    }

    const StringsList allowed_vel_init_methods = { "random", "manual" };

    if (!StringUtils::labelInArray(init_vel_type, allowed_vel_init_methods))
        throw std::invalid_argument(std::format("The specified velocity initialization method ({}) is not supported!",
            init_vel_type));

    if (init_vel_type == "manual" && !init_vel_specification.empty()) {
        try {
            const int dummy = 0;
            // Try using specification as filename format
            std::string formatted_filename = std::vformat(init_vel_specification, std::make_format_args(dummy));
            /// TODO: What if both 0 and 1 don't exist?
            config.init_vel_index_offset = static_cast<int>(!std::filesystem::exists(formatted_filename));
            config.init_vel_filename = init_vel_specification;
        } catch (const std::format_error&) {
            throw std::invalid_argument(
                std::format("The filename format ({}) for velocity initialization is invalid!",
                    init_vel_specification)
            );
        } catch (...) {
            throw std::runtime_error(
                std::format("Filename format ({}) for velocity initialization validation failed",
                    init_vel_specification)
            );
        }
    }

    /// TODO: Local init_vel_type isn't really necessary
    config.init_vel_type = init_vel_type;
}

/**
 * Load parameters pertaining to the external potential.
 *
 * @param config Config object to load parameters into.
 */
void Params::loadExternalPotentialParams(SimulationConfig& config) const {
    const StringsList allowed_ext_potential_names = {"free", "harmonic", "double_well", "cosine"};

    config.ext_pot_name = m_reader.GetString(Sections::EXT_POTENTIAL, "name", "free");
    if (!StringUtils::labelInArray(config.ext_pot_name, allowed_ext_potential_names))
        throw std::invalid_argument(std::format("The specified external potential ({}) is not supported!",
            config.ext_pot_name));

    /// TODO: Convert to switch-case
    if (config.ext_pot_name == "harmonic") {
        // In atomic units, the angular frequency of the oscillator has the same dimensions as the energy
        config.ext_pot_params["omega"] = Units::getQuantity(
            "energy", m_reader.Get(Sections::EXT_POTENTIAL, "omega", "1.0 millielectronvolt"));
    } else if (config.ext_pot_name == "double_well") {
        config.ext_pot_params["strength"] = Units::getQuantity(
            "energy", m_reader.Get(Sections::EXT_POTENTIAL, "strength", "1.0 millielectronvolt"));
        config.ext_pot_params["location"] = Units::getQuantity(
            "length", m_reader.Get(Sections::EXT_POTENTIAL, "location", "1.0 angstrom"));
    } else if (config.ext_pot_name == "cosine") {
        config.ext_pot_params["amplitude"] = Units::getQuantity(
            "energy", m_reader.Get(Sections::EXT_POTENTIAL, "amplitude", "1.0 millielectronvolt"));
        config.ext_pot_params["phase"] = m_reader.GetReal(Sections::EXT_POTENTIAL, "phase", 1.0);
    }
}

/**
 * Load parameters pertaining to the interaction potential.
 *
 * @param config Config object to load parameters into.
 */
void Params::loadInteractionPotentialParams(SimulationConfig& config) const {
    const StringsList allowed_int_potential_names = { "aziz", "free", "harmonic", "dipole" };

    config.int_pot_name = m_reader.GetString(Sections::INT_POTENTIAL, "name", "free");

    if (!StringUtils::labelInArray(config.int_pot_name, allowed_int_potential_names))
        throw std::invalid_argument(std::format("The specified interaction potential ({}) is not supported!",
            config.int_pot_name));

    config.int_pot_params["cutoff"] = Units::getQuantity("length", m_reader.Get(Sections::INT_POTENTIAL, "cutoff", "-1.0 angstrom"));

    /// TODO: Convert to switch-case
    if (config.int_pot_name == "free") {
        // In the special case of free particles, the cutoff distance is set to zero
        config.int_pot_params["cutoff"] = 0.0;
    } else if (config.int_pot_name == "harmonic") {
        // In atomic units, the angular frequency of the oscillator has the same dimensions as the energy
        config.int_pot_params["omega"] = Units::getQuantity(
            "energy", m_reader.Get(Sections::INT_POTENTIAL, "omega", "1.0 millielectronvolt"));
    } else if (config.int_pot_name == "double_well") {
        config.int_pot_params["strength"] = Units::getQuantity(
            "energy", m_reader.Get(Sections::INT_POTENTIAL, "strength", "1.0 millielectronvolt"));
        config.int_pot_params["location"] = Units::getQuantity(
            "length", m_reader.Get(Sections::INT_POTENTIAL, "location", "1.0 angstrom"));
    } else if (config.int_pot_name == "dipole") {
        config.int_pot_params["strength"] = m_reader.GetReal(Sections::INT_POTENTIAL, "strength", 1.0);
    }
}

/**
 * Load parameters pertaining to the output (dumps).
 *
 * @param config Config object to load parameters into.
 */
void Params::loadOutputParams(SimulationConfig& config) const {
    config.dumps_list["positions"] = m_reader.Get(Sections::DUMP, "positions", "off");
    config.dumps_list["velocities"] = m_reader.Get(Sections::DUMP, "velocities", "off");
    config.dumps_list["forces"] = m_reader.Get(Sections::DUMP, "forces", "off");
}

/**
 * Load parameters pertaining to the observables.
 *
 * @param config Config object to load parameters into.
 */
void Params::loadObservableParams(SimulationConfig& config) const {
    config.observables_list["energy"] = m_reader.Get(Sections::OBSERVABLES, "energy", "kelvin");
    config.observables_list["classical"] = m_reader.Get(Sections::OBSERVABLES, "classical", "off");
    config.observables_list["bosonic"] = m_reader.Get(Sections::OBSERVABLES, "bosonic", "off");
    config.observables_list["gsf"] = m_reader.Get(Sections::OBSERVABLES, "gsf", "off");
}
