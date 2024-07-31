#include "params.h"

#include <sstream>
#include <regex>
#include <format>

Params::Params(const std::string& filename, const int& rank) : reader(filename) {
    printStatus("Initializing the simulation parameters", rank);

    if (reader.ParseError() < 0)
        throw std::invalid_argument(std::format("Unable to read the configuration file {}", filename));

    /****** Simulation params ******/
    sim["dt"] = getQuantity("time", reader.Get(Sections::SIMULATION, "dt", "1.0 femtosecond"));
    sim["threshold"] = reader.GetReal(Sections::SIMULATION, "threshold", 0.1);
    sim["gamma"] = reader.GetReal(Sections::SIMULATION, "gamma", -1.0);

    if (std::get<double>(sim["gamma"]) < 0)
        sim["gamma"] = 1 / (100.0 * std::get<double>(sim["dt"]));

    sim["steps"] = static_cast<long>(
        std::stod(reader.Get(Sections::SIMULATION, "steps", "1e5")));  // Scientific notation
    sim["sfreq"] = reader.GetLong(Sections::SIMULATION, "sfreq", 1000); /// @todo: Add support for scientific notation
    sim["enable_t"] = reader.GetBoolean(Sections::SIMULATION, "enable_thermostat", true);
    sim["nbeads"] = reader.GetInteger(Sections::SIMULATION, "nbeads", 4);

    // if (int nbeads = std::get<int>(sim["nbeads"]); nbeads < 2)
    //     throw std::invalid_argument(std::format("The specified number of beads ({}) is less than two!", nbeads));

    if (int nbeads = std::get<int>(sim["nbeads"]); nbeads < 1)
        throw std::invalid_argument(std::format("The specified number of beads ({}) is less than one!", nbeads));

    sim["seed"] = static_cast<unsigned int>(std::stod(reader.Get(Sections::SIMULATION, "seed", "1234")));

    /****** Flags ******/
    // Bosonic or distinguishable simulation?
    sim["bosonic"] = reader.GetBoolean(Sections::SIMULATION, "bosonic", false);
    // Fix the center of mass?
    sim["fixcom"] = reader.GetBoolean(Sections::SIMULATION, "fixcom", true);
    // Enable periodic boundary conditions?
    sim["pbc"] = reader.GetBoolean(Sections::SIMULATION, "pbc", false);

    sim["apply_mic_spring"] = reader.GetBoolean(Sections::SIMULATION, "apply_mic_spring", false);
    sim["apply_mic_potential"] = reader.GetBoolean(Sections::SIMULATION, "apply_mic_potential", true);
    sim["apply_wrap"] = reader.GetBoolean(Sections::SIMULATION, "apply_wrap", false);
    sim["apply_wrap_first"] = reader.GetBoolean(Sections::SIMULATION, "apply_wrap_first", false);
    sim["apply_wind"] = reader.GetBoolean(Sections::SIMULATION, "apply_wind", false);

    // Maximum winding sector to consider for periodic boundary conditions
    sim["max_wind"] = reader.GetInteger(Sections::SIMULATION, "max_wind", 0);

    std::string init_pos_type, init_pos_specification;

    if (!parseTokenParentheses(reader.Get(Sections::SIMULATION, "initial_position", "random"), init_pos_type,
                               init_pos_specification)) {
        throw std::invalid_argument("The coordinate initialization method format is invalid!");
    }

    allowed_coord_init_methods = { "random", "xyz", "grid" }; /// @todo Add "cell" option

    if (!labelInArray(init_pos_type, allowed_coord_init_methods))
        throw std::invalid_argument(std::format("The specified coordinate initialization method ({}) is not supported!",
                                                init_pos_type));

    sim["init_pos_type"] = init_pos_type;

    if (init_pos_type == "xyz") {
        sim["init_pos_xyz_filename"] = init_pos_specification;
    }

    // Allowed velocity initialization methods:
    // "random": samples from Maxwell-Boltzmann distribution
    // "manual": reads from vel_X.dat files
    allowed_vel_init_methods = { "random", "manual" };
    std::string init_vel_type = reader.GetString(Sections::SIMULATION, "initial_velocity", "random");
    sim["init_vel_type"] = init_vel_type;

    if (!labelInArray(init_vel_type, allowed_vel_init_methods))
        throw std::invalid_argument(std::format("The specified velocity initialization method ({}) is not supported!",
                                                init_vel_type));

    /* System params */
    sys["temperature"] = getQuantity("temperature", reader.Get(Sections::SYSTEM, "temperature", "1.0 kelvin"));
    if (double temp = std::get<double>(sys["temperature"]); temp <= 0.0) {
        throw std::invalid_argument(std::format("The specified temperature ({0:4.3f} kelvin) is unphysical!", temp));
    }

    sys["natoms"] = reader.GetInteger(Sections::SYSTEM, "natoms", 1);
    if (int natoms = std::get<int>(sys["natoms"]); natoms < 1)
        throw std::invalid_argument(std::format("The specified number of particles ({}) is smaller than one!", natoms));

    sys["mass"] = getQuantity("mass", reader.Get(Sections::SYSTEM, "mass", "1.0 dalton"));
    if (double mass = std::get<double>(sys["mass"]); mass <= 0.0)
        throw std::invalid_argument(std::format("The provided mass ({0:4.3f}) is unphysical!", mass));

    sys["size"] = getQuantity("length", reader.Get(Sections::SYSTEM, "size", "1.0 picometer"));
    if (double size = std::get<double>(sys["size"]); size <= 0.0)
        throw std::invalid_argument(std::format("The provided system size ({0:4.3f}) is unphysical!", size));

    allowed_int_potential_names = { "aziz", "free", "harmonic", "dipole" };
    allowed_ext_potential_names = { "free", "harmonic", "double_well" };

    /****** Interaction potential ******/

    std::string interaction_name = reader.GetString(Sections::INT_POTENTIAL, "name", "free");

    if (!labelInArray(interaction_name, allowed_int_potential_names))
        throw std::invalid_argument(std::format("The specified interaction potential ({}) is not supported!",
                                                interaction_name));

    interaction_pot["name"] = interaction_name;
    interaction_pot["cutoff"] = getQuantity("length", reader.Get(Sections::INT_POTENTIAL, "cutoff", "-1.0 angstrom"));


    if (interaction_name == "free") {
        // In the special case of free particles, the cutoff distance is set to zero
        interaction_pot["cutoff"] = 0.0;
    } else if (interaction_name == "harmonic") {
        // In atomic units, the angular frequency of the oscillator has the same dimensions as the energy
        interaction_pot["omega"] = getQuantity(
            "energy", reader.Get(Sections::INT_POTENTIAL, "omega", "1.0 millielectronvolt"));
    } else if (interaction_name == "double_well") {
        interaction_pot["strength"] = getQuantity(
            "energy", reader.Get(Sections::INT_POTENTIAL, "strength", "1.0 millielectronvolt"));
        interaction_pot["location"] = getQuantity(
            "length", reader.Get(Sections::INT_POTENTIAL, "location", "1.0 angstrom"));
    } else if (interaction_name == "dipole") {
        interaction_pot["strength"] = reader.GetReal(Sections::INT_POTENTIAL, "strength", 1.0);
    }

    /****** External potential ******/

    std::string external_name = reader.GetString(Sections::EXT_POTENTIAL, "name", "free");

    if (!labelInArray(external_name, allowed_ext_potential_names))
        throw std::invalid_argument(std::format("The specified external potential ({}) is not supported!",
                                                external_name));

    external_pot["name"] = external_name;

    if (external_name == "harmonic") {
        // In atomic units, the angular frequency of the oscillator has the same dimensions as the energy
        external_pot["omega"] = getQuantity(
            "energy", reader.Get(Sections::EXT_POTENTIAL, "omega", "1.0 millielectronvolt"));
    } else if (external_name == "double_well") {
        external_pot["strength"] = getQuantity(
            "energy", reader.Get(Sections::EXT_POTENTIAL, "strength", "1.0 millielectronvolt"));
        external_pot["location"] = getQuantity(
            "length", reader.Get(Sections::EXT_POTENTIAL, "location", "1.0 angstrom"));
    }

    /****** Output ******/

    out["positions"] = reader.GetBoolean(Sections::OUTPUT, "positions", false);
    out["velocities"] = reader.GetBoolean(Sections::OUTPUT, "velocities", false);
    out["forces"] = reader.GetBoolean(Sections::OUTPUT, "forces", false);
    out["wind_prob"] = reader.GetBoolean(Sections::OUTPUT, "wind_prob", false);
}

bool Params::labelInArray(const std::string& label, const StringsList& arr) {
    return std::ranges::find(arr, label) != arr.end();
}

// Splits a string into several strings based on a delimiter
std::vector<std::string> Params::splitString(const std::string& line, char delimiter) {
    std::stringstream ss(line);
    std::vector<std::string> tokens;
    std::string token;

    while (std::getline(ss, token, delimiter)) {
        if (!token.empty())
            tokens.push_back(token);
    }

    return tokens;
}

// Parse a string containing a numerical value and a unit
std::pair<double, std::string> Params::parseQuantity(const std::string& input) {
    std::istringstream iss(input);
    std::string unit;

    // Read the floating-point number and the unit
    // See: https://en.cppreference.com/w/cpp/io/basic_istream/operator_gtgt
    if (double value; iss >> value >> std::ws >> unit) {
        return std::make_pair(value, unit);
    }

    throw std::invalid_argument("Invalid input format");
}

/**
 * Accepts a string of the format "<number> <unit>" and returns the numerical
 * value of the quantity in internal (atomic) units.
 *
 * @param family The family of units to which the quantity belongs.
 * @param input The string containing the numerical value and the unit.
 * @return The numerical value of the quantity in internal units.
 */
double Params::getQuantity(const std::string& family, const std::string& input) {
    // Extract numerical value and specified unit
    auto [value, raw_unit] = parseQuantity(input);

    // Match the provided unit to known units.
    return Units::convertToInternal(family, raw_unit, value);
}

// Used to parse strings of the format "token(value)"
bool Params::parseTokenParentheses(const std::string& input, std::string& token, std::string& value) {
    // Define a regular expression for the specified format
    // @todo Generalize regex to catch cases such as "foo()" and "foo" (without parentheses)	
    // Catch "foo(bar)", "foo()" and "foo"
    const std::regex pattern(R"(\s*(\w+)\((.*)\)\s*|\s*(\w+)\s*)");

    // Match the input string against the pattern
    if (std::smatch matches; std::regex_match(input, matches, pattern)) {
        // Check if there are matched groups
        if (matches.size() == 4 && (matches[1].matched || matches[3].matched)) {
            // Extract token and value from the matched groups
            token = (matches[1].matched) ? matches[1].str() : matches[3].str();
            value = matches[2].str();
            return true;
        }
    }

    // If no match is found or the match does not have the expected groups, return false
    return false;
}
