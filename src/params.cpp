#include "params.h"
#include "core/potential_config.h"
#include "../libs/string_utils.h"

#include <regex>
#include <format>
#include <filesystem>

using Quantity = Units::Quantity;

namespace {
    /**
     * Shared validation logic for "xyz" / "xyz(format)" style position/velocity initialization.
     * Verifies the (optional) filename format, determines whether on-disk bead numbering starts at
     * 0 or 1, and enforces the mandatory unit and frame-selection parameters.
     *
     * @param reader          The INI reader to pull the unit/frame parameters from.
     * @param section         Section to read the unit/frame parameters from.
     * @param label           Human-readable label used in error messages (e.g. "Coordinate", "Velocity").
     * @param key_prefix      Key prefix shared by the unit/frame parameters (e.g. "initial_position").
     * @param specification   The filename/format string parsed out of "xyz(...)" syntax (may be empty).
     * @param[out] index_offset Set to 1 if bead index 0 is missing on disk (numbering starts at 1), else 0.
     * @param[out] filename     Set to the (unformatted) filename/format string, stored verbatim.
     * @param[out] unit         Set to the user-specified unit string.
     * @param[out] frame        Set to the selected frame index.
     * @param[out] frame_mode   Set to the selected frame-selection mode.
     */
    template <typename Reader>
    void loadXyzInitParams(
        const Reader& reader,
        const std::string& section,
        const std::string& label,
        const std::string& key_prefix,
        const std::string& specification,
        int& index_offset,
        std::string& filename,
        std::string& unit,
        long& frame,
        XyzFrameSelectionMode& frame_mode
    ) {
        try {
            // If the user provided a filename format, try to format it using the first bead index
            constexpr int dummy = 0;
            const std::string formatted_filename = std::vformat(specification, std::make_format_args(dummy));
            /// TODO: What if both 0 and 1 don't exist?
            index_offset = static_cast<int>(!std::filesystem::exists(formatted_filename));
            filename = specification;
        } catch (const std::format_error&) {
            throw std::invalid_argument(
                std::format("The filename format ({}) for {} initialization is invalid!", specification, label));
        } catch (...) {
            throw std::runtime_error(
                std::format("Filename format ({}) for {} initialization validation failed", specification, label));
        }

        // If the initialization method is "xyz" but the user didn't provide a unit,
        // we throw an error instead of using a default (don't be smarter than the user)
        const std::string unit_key = key_prefix + "_unit";
        if (!reader.HasValue(section, unit_key)) {
            throw std::invalid_argument(
                label + " initialization method is 'xyz' but the units have not been specified (use '" + unit_key + "')");
        }
        unit = reader.Get(section, unit_key, "");

        frame = reader.GetLong(section, key_prefix + "_frame", 0);
        if (frame < 0) {
            throw std::invalid_argument(key_prefix + "_frame must be non-negative");
        }

        const std::string frame_mode_str = reader.GetString(section, key_prefix + "_frame_mode", "index");
        if (frame_mode_str == "index") {
            frame_mode = XyzFrameSelectionMode::Index;
        } else if (frame_mode_str == "step") {
            frame_mode = XyzFrameSelectionMode::Step;
        } else {
            throw std::invalid_argument(std::format(
                "Unsupported {}_frame_mode '{}'; expected 'index' or 'step'", key_prefix, frame_mode_str));
        }
    }

}

Params::Params(const std::string& filename, int rank) : m_reader(filename), m_rank(rank) {
    printStatus("Initializing the simulation parameters", rank);

    if (m_reader.ParseError() < 0)
        throw std::invalid_argument(std::format("Unable to read the configuration file {}", filename));
}

std::shared_ptr<SimulationConfig> Params::load() const {
    auto config = std::make_shared<SimulationConfig>();

    config->this_bead = m_rank;

    // TODO: Check for unknown headers. Warn the user if any are found.


    // NOTE: The calls below have implicit ordering dependencies - reordering them can silently
    // leave config fields unset or unvalidated. Known dependencies:
    //   - loadPropagatorParams needs config.bosonic (set by loadSimulationParams)
    //   - loadThermostatParams needs config.mass (set by loadSystemParams) and config.dt,
    //     config.nbeads (set by loadSimulationParams)
    //   - loadExternalPotentialParams / loadInteractionPotentialParams need config.mass and
    //     config.box_size (set by loadSystemParams)
    // If you add or reorder a load*Params function, please keep this list in sync.
    loadSimulationParams(*config);
    loadSystemParams(*config);
    loadPropagatorParams(*config);
    loadThermostatParams(*config);
    loadCoordInitParams(*config);
    loadVelocityInitParams(*config);
    loadExternalPotentialParams(*config);
    loadInteractionPotentialParams(*config);
    loadOutputParams(*config);
    loadObservableParams(*config);
    loadRpmdParams(*config);

    return config;
}

double Params::loadQuantity(
    const std::string& family,
    const std::string& section,
    const std::string& key,
    const std::string& default_value,
    double& destination,
    StringMap& units_destination
) const {
    const Quantity quantity = Quantity::create(
        family,
        m_reader.Get(section, key, default_value)
    );

    destination = quantity.toInternal();
    units_destination[key] = quantity.unit;

    return quantity.value; // Return the original user-provided value (not converted to internal units)
}

/**
 * Load basic parameters pertaining to the simulation.
 *
 * @param config Config object to load parameters into.
 */
void Params::loadSimulationParams(SimulationConfig& config) const {
    //config.dt = Units::getQuantity("time", m_reader.Get(Sections::SIMULATION, "dt", "1.0 femtosecond"));
    const double dt = loadQuantity(
        "time", 
        Sections::SIMULATION, 
        "dt", 
        "1.0 femtosecond", 
        config.dt, 
        config.units_list
    );

    if (dt <= 0.0)
        throw std::invalid_argument(std::format("Invalid time-step ({} is not positive)", dt));

    // Config's threshold is specified as a fraction of the total number of steps,
    // so we shall convert it to an absolute number of steps later on.
    // For now, we just validate that the user-provided threshold is in [0, 1].
    const double threshold = m_reader.GetReal(Sections::SIMULATION, "threshold", 0.0);
    if (threshold < 0.0 || threshold > 1.0)
        throw std::invalid_argument(std::format("Invalid threshold ({} must lie in [0, 1])", threshold));

    config.steps = static_cast<long>(std::stod(m_reader.Get(Sections::SIMULATION, "steps", "1e5")));
    if (config.steps < 1)
        throw std::invalid_argument(std::format("Invalid number of steps ({}<1)", config.steps));

    config.threshold = static_cast<long>(threshold * config.steps);  // Threshold in terms of steps

    //config.sfreq = m_reader.GetLong(Sections::SIMULATION, "sfreq", 1000); /// @todo: Add support for scientific notation
    config.sfreq = static_cast<long>(std::stod(m_reader.Get(Sections::SIMULATION, "sfreq", "1e3")));
    if (config.sfreq < 1)
        throw std::invalid_argument(std::format("Invalid value for observable recording frequency ({}<1)", config.sfreq));

    config.nbeads = m_reader.GetInteger(Sections::SIMULATION, "nbeads", 4);
    if (config.nbeads < 1)
        throw std::invalid_argument(std::format("Invalid number of beads ({}<1)", config.nbeads));

    // unsigned int seed = static_cast<unsigned int>(time(nullptr))
    // TODO: For large seed values, we might be at risk of a silent precision loss through the double round-trip;
    // GetInteger/GetLong-style integer parsing (or a small helper that accepts scientific notation without going through double) would be safer, and
    // would match the pattern already used for steps/sfreq.
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

    //config.mass = Units::getQuantity("mass", m_reader.Get(Sections::SYSTEM, "mass", "1.0 dalton"));
    const double mass = loadQuantity(
        "mass",
        Sections::SYSTEM,
        "mass",
        "1.0 dalton",
        config.mass,
        config.units_list
    );

    if (mass <= 0.0)
        throw std::invalid_argument(std::format("The provided mass ({0:4.3f}) is unphysical!", mass));

    //config.box_size = Units::getQuantity("length", m_reader.Get(Sections::SYSTEM, "size", "1.0 picometer"));
    const double box_size = loadQuantity(
        "length",
        Sections::SYSTEM,
        "size",
        "1.0 picometer",
        config.box_size,
        config.units_list
    );

    if (box_size <= 0.0)
        throw std::invalid_argument(std::format("The provided system size ({0:4.3f}) is unphysical!", box_size));
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
    using namespace std::string_view_literals;
    constexpr auto allowed_propagators = std::array{ "cartesian"sv, "normal_modes"sv };
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
    using namespace std::string_view_literals;
    constexpr auto allowed_thermostats = std::array{
        "langevin"sv,            // A Langevin thermostat coupled to the Cartesian coordinates
        "nose_hoover"sv,         // A single Nose-Hoover chain coupled to the whole system
        "nose_hoover_np"sv,      // A unique Nose-Hoover chain coupled to each particle
        "nose_hoover_np_dim"sv,  // A unique Nose-Hoover chain coupled to each Cartesian coordinate of each particle
        "none"sv                 // No thermostat (NVE simulation)
    };

    /// TODO: Perhaps move config.rpmd_config.enabled initialization to here.
    ///       Pros: We can allow the user to not specify the thermostat type
    ///       when RPMD is enabled (set it to "none" automatically.)
    ///       Cons: Moves RPMD-specific logic into a thermostat-specific function, which may not be ideal.
    ///       Also, it may be better to require the user to explicitly specify the thermostat type even when RPMD is
    ///       enabled, for clarity.

    // Read and validate thermostat type first (this drives other validations)
    config.thermostat.type = m_reader.GetString(Sections::SIMULATION, "thermostat", "error");
    if (config.thermostat.type == "error") {
        throw std::invalid_argument("Thermostat must be specified!");
    }
    if (!StringUtils::labelInArray(config.thermostat.type, allowed_thermostats)) {
        throw std::invalid_argument(std::format("Unsupported thermostat type ({})", config.thermostat.type));
    }

    // Validate temperature. When no thermostat is used, temperature value may still be used for initializing velocities
    //config.temperature = Units::getQuantity("temperature", m_reader.Get(Sections::SYSTEM, "temperature", "1.0 kelvin"));
    const double temperature = loadQuantity(
        "temperature",
        Sections::SYSTEM,
        "temperature",
        "1.0 kelvin",
        config.temperature,
        config.units_list
    );

    if (temperature <= 0) {
        throw std::invalid_argument(
            std::format(
                "The specified temperature ({0:4.3f} {1}) is not positive", 
                temperature,
                config.units_list["temperature"]
            )
        );
    }

    // Define additional parameters based on rudimentary config parameters
    config.beta = 1.0 / (Constants::kB * config.temperature);

#if TAU_CONVENTION
    // Tau convention [J. Chem. Phys. 133, 124104 (2010)]
    config.thermo_beta = config.beta / config.nbeads;
    config.omega_p = config.nbeads / (config.beta * Constants::hbar);
#else
    // Beta convention
    config.thermo_beta = config.beta;
    config.omega_p = sqrt(config.nbeads) / (config.beta * Constants::hbar);
#endif

    config.spring_constant = config.mass * config.omega_p * config.omega_p;
    config.beta_half_k = config.thermo_beta * 0.5 * config.spring_constant;

    // Determine if this is a Nose-Hoover type thermostat
    const bool is_nose_hoover = config.thermostat.type.find("nose_hoover") != std::string::npos;
    config.thermostat.is_nose_hoover = is_nose_hoover;

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
        config.thermostat.nchains = nchains;
    }

    // Handle normal mode thermostat coupling
    bool nmthermostat = m_reader.GetBoolean(Sections::SIMULATION, "nmthermostat", false);
    if (nmthermostat && config.thermostat.type == "none") {
        throw std::invalid_argument("nmthermostat cannot be used in NVE ensemble!");
    }

    // Only set nmthermostat if we're using a thermostat
    config.thermostat.couple_to_nm = nmthermostat;

    // Handle gamma parameter (only relevant for Langevin thermostat)
    if (config.thermostat.type == "langevin") {
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

        config.thermostat.gamma = gamma;
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

    using namespace std::string_view_literals;
    constexpr auto allowed_coord_init_methods = std::array{ "random"sv, "xyz"sv, "grid"sv }; /// TODO: Use cell instead of grid?

    if (!StringUtils::labelInArray(init_pos_type, allowed_coord_init_methods))
        throw std::invalid_argument(std::format("The specified coordinate initialization method ({}) is not supported!",
            init_pos_type));

    if (init_pos_type == "xyz") {
        /// TODO: Check correctness of the unit (perhaps create a universal function for this, and do this in other places as well)
        loadXyzInitParams(
            m_reader,
            Sections::SIMULATION,
            "Coordinate",
            "initial_position",
            init_pos_specification,
            config.init_pos_index_offset,
            config.init_pos_filename,
            config.init_pos_unit,
            config.init_pos_frame,
            config.init_pos_frame_mode
        );
    } else if (
        m_reader.HasValue(Sections::SIMULATION, "initial_position_frame") ||
        m_reader.HasValue(Sections::SIMULATION, "initial_position_frame_mode")
    ) {
        throw std::invalid_argument(
            "initial_position_frame and initial_position_frame_mode can only be used with initial_position = xyz(...)"
        );
    }

    config.init_pos_type = init_pos_type;
}

/**
 * Load parameters pertaining to the velocity initialization.
 *
 * @param config Config object to load parameters into.
 */
void Params::loadVelocityInitParams(SimulationConfig& config) const {
    // Allowed velocity initialization methods:
    //  "random": samples from Maxwell-Boltzmann distribution
    //  "xyz": reads from strict XYZ format files
    //  "xyz(format)": reads from XYZ format with per-bead files
    std::string init_vel_type, init_vel_specification;

    if (!StringUtils::parseTokenParentheses(
        m_reader.Get(Sections::SIMULATION, "initial_velocity", "random"), 
        init_vel_type,
        init_vel_specification
        )) {
        throw std::invalid_argument("Invalid velocity initialization method");
    }

    using namespace std::string_view_literals;
    constexpr auto allowed_vel_init_methods = std::array{ "random"sv, "xyz"sv };

    if (!StringUtils::labelInArray(init_vel_type, allowed_vel_init_methods))
        throw std::invalid_argument(std::format("The specified velocity initialization method ({}) is not supported!",
            init_vel_type));

    if (init_vel_type == "xyz") {
        loadXyzInitParams(
            m_reader,
            Sections::SIMULATION,
            "Velocity",
            "initial_velocity",
            init_vel_specification,
            config.init_vel_index_offset,
            config.init_vel_filename,
            config.init_vel_unit,
            config.init_vel_frame,
            config.init_vel_frame_mode
        );
    } else if (
        m_reader.HasValue(Sections::SIMULATION, "initial_velocity_frame") ||
        m_reader.HasValue(Sections::SIMULATION, "initial_velocity_frame_mode")) {
        throw std::invalid_argument(
            "initial_velocity_frame and initial_velocity_frame_mode can only be used with initial_velocity = xyz(...)"
        );
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
    using namespace std::string_view_literals;
    constexpr auto allowed_ext_potential_names = std::array{"free"sv, "harmonic"sv, "anharmonic"sv, "double_well"sv, "cosine"sv};

    const std::string name = m_reader.GetString(Sections::EXT_POTENTIAL, "name", "free");
    if (!StringUtils::labelInArray(name, allowed_ext_potential_names))
        throw std::invalid_argument(std::format("The specified external potential ({}) is not supported!", name));

    if (name == "harmonic") {
        // In atomic units, the angular frequency of the oscillator has the same dimensions as the energy
        const double omega = Units::getQuantity(
            "energy", m_reader.Get(Sections::EXT_POTENTIAL, "omega", "1.0 millielectronvolt"));
        config.ext_potential_cfg = PotentialConfig{name, HarmonicPotentialParams{config.mass, omega}};
    } else if (name == "anharmonic") {
        const double omega = Units::getQuantity(
            "energy", m_reader.Get(Sections::EXT_POTENTIAL, "omega", "1.0 millielectronvolt"));
        const double cubic_const = m_reader.GetReal(Sections::EXT_POTENTIAL, "cubic_const", 0.0);
        const double quart_const = m_reader.GetReal(Sections::EXT_POTENTIAL, "quart_const", 0.0);
        config.ext_potential_cfg = PotentialConfig{name, AnharmonicPotentialParams{config.mass, omega, cubic_const, quart_const}};
    } else if (name == "double_well") {
        const double strength = Units::getQuantity(
            "energy", m_reader.Get(Sections::EXT_POTENTIAL, "strength", "1.0 millielectronvolt"));
        const double location = Units::getQuantity(
            "length", m_reader.Get(Sections::EXT_POTENTIAL, "location", "1.0 angstrom"));
        config.ext_potential_cfg = PotentialConfig{name, DoubleWellPotentialParams{config.mass, strength, location}};
    } else if (name == "cosine") {
        const double amplitude = Units::getQuantity(
            "energy", m_reader.Get(Sections::EXT_POTENTIAL, "amplitude", "1.0 millielectronvolt"));
        const double phase = m_reader.GetReal(Sections::EXT_POTENTIAL, "phase", 0.0);
        config.ext_potential_cfg = PotentialConfig{name, CosinePotentialParams{amplitude, config.box_size, phase}};
    } else {
        config.ext_potential_cfg = PotentialConfig{};  // defaults to "free"
    }
}

/**
 * Load parameters pertaining to the interaction potential.
 *
 * @param config Config object to load parameters into.
 */
void Params::loadInteractionPotentialParams(SimulationConfig& config) const {
    using namespace std::string_view_literals;
    constexpr auto allowed_int_potential_names = std::array{ "aziz"sv, "free"sv, "harmonic"sv, "anharmonic"sv, "dipole"sv };

    const std::string name = m_reader.GetString(Sections::INT_POTENTIAL, "name", "free");
    if (!StringUtils::labelInArray(name, allowed_int_potential_names))
        throw std::invalid_argument(std::format("The specified interaction potential ({}) is not supported!", name));

    // In the special case of free particles the cutoff is forced to zero (no pairwise loops).
    // For all other potentials the user-supplied cutoff is used (negative = no cutoff).
    const double cutoff = (name == "free")
        ? 0.0
        : Units::getQuantity("length", m_reader.Get(Sections::INT_POTENTIAL, "cutoff", "-1.0 angstrom"));

    if (name == "free") {
        config.int_potential_cfg = PotentialConfig{name, FreePotentialParams{}, cutoff};
    } else if (name == "harmonic") {
        const double omega = Units::getQuantity(
            "energy", m_reader.Get(Sections::INT_POTENTIAL, "omega", "1.0 millielectronvolt"));
        config.int_potential_cfg = PotentialConfig{name, HarmonicPotentialParams{config.mass, omega}, cutoff};
    } else if (name == "anharmonic") {
        const double omega = Units::getQuantity(
            "energy", m_reader.Get(Sections::INT_POTENTIAL, "omega", "1.0 millielectronvolt"));
        const double cubic_const = m_reader.GetReal(Sections::INT_POTENTIAL, "cubic_const", 0.0);
        const double quart_const = m_reader.GetReal(Sections::INT_POTENTIAL, "quart_const", 0.0);
        config.int_potential_cfg = PotentialConfig{name, AnharmonicPotentialParams{config.mass, omega, cubic_const, quart_const}, cutoff};
    } else if (name == "double_well") {
        const double strength = Units::getQuantity(
            "energy", m_reader.Get(Sections::INT_POTENTIAL, "strength", "1.0 millielectronvolt"));
        const double location = Units::getQuantity(
            "length", m_reader.Get(Sections::INT_POTENTIAL, "location", "1.0 angstrom"));
        config.int_potential_cfg = PotentialConfig{name, DoubleWellPotentialParams{config.mass, strength, location}, cutoff};
    } else if (name == "dipole") {
        const double strength = m_reader.GetReal(Sections::INT_POTENTIAL, "strength", 1.0);
        config.int_potential_cfg = PotentialConfig{name, DipolePotentialParams{strength}, cutoff};
    } else if (name == "aziz") {
        config.int_potential_cfg = PotentialConfig{name, AzizPotentialParams{}, cutoff};
    } else {
        config.int_potential_cfg = PotentialConfig{name, FreePotentialParams{}, cutoff};
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
    config.observables_list["center_of_mass"] = m_reader.Get(Sections::OBSERVABLES, "center_of_mass", "off");
}


void Params::loadRpmdParams(SimulationConfig& config) const {
    config.rpmd_config.enabled = m_reader.GetBoolean(Sections::SIMULATION, "rpmd", false);

    if (config.rpmd_config.enabled) {
        if (config.thermostat.type != "none") {
            throw std::invalid_argument(
                std::format(
                    "RPMD mode requires 'thermostat = none' but got '{}'", 
                    config.thermostat.type
                )
            );
        }

        if (
            m_reader.HasValue(Sections::SIMULATION, "initial_position_frame") ||
            m_reader.HasValue(Sections::SIMULATION, "initial_position_frame_mode") ||
            m_reader.HasValue(Sections::SIMULATION, "initial_velocity_frame") ||
            m_reader.HasValue(Sections::SIMULATION, "initial_velocity_frame_mode")
            )
        {
            throw std::invalid_argument(
                "RPMD cannot be used simultaneously with 'initial_position_frame', 'initial_position_frame_mode', "
                "'initial_velocity_frame', or 'initial_velocity_frame_mode' parameters."
            );
        }

        config.rpmd_config.num_runs = m_reader.GetInteger(Sections::SIMULATION, "rpmd_nruns", 1);
        if (config.rpmd_config.num_runs < 1)
            throw std::invalid_argument("rpmd_nruns must be at least 1");

        config.rpmd_config.nvt_discard_frac = m_reader.GetReal(Sections::SIMULATION, "rpmd_nvt_discard_fraction", 0.0);
        if (config.rpmd_config.nvt_discard_frac < 0 || config.rpmd_config.nvt_discard_frac >= 1)
            throw std::invalid_argument("rpmd_nvt_discard_fraction must be non-negative and less than 1");

        /*config.rpmd_config.xyz_path = m_reader.Get(Sections::SIMULATION, "rpmd_xyz_path", "");
        if (config.rpmd_config.xyz_path.empty())
            throw std::invalid_argument("rpmd_xyz_path must be specified when RPMD is enabled");*/
    }
}
