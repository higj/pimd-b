#pragma once

#include "common.h"
#include "inireader.h"

#include <vector>
#include <unordered_map>
#include <variant>

namespace Sections {
    const std::string SYSTEM = "system";
    const std::string SIMULATION = "simulation";
    const std::string EXT_POTENTIAL = "external_potential";
    const std::string INT_POTENTIAL = "interaction_potential";
    const std::string OUTPUT = "output";
    const std::string OBSERVABLES = "observables";
}

class Params {
public:
    explicit Params(const std::string& filename, const int& rank);
    INIReader reader;

    // Map holding the simulation parameters
    VariantMap sim;
    // Map holding the (physical) system parameters
    VariantMap sys;
    // Map holding the interaction potential parameters
    VariantMap interaction_pot;
    // Map holding the interaction potential parameters
    VariantMap external_pot;
    // Map holding the output parameters
    StringMap states;
    // Map holding the observable settings
    StringMap observables;

    /// @todo Consider refactoring this method to a more general utility class (maybe to units)?
    static std::pair<double, std::string> parseQuantity(const std::string& input);
    static double getQuantity(const std::string& family, const std::string& input);
    static bool parseTokenParentheses(const std::string& input, std::string& token, std::string& value);

private:
    StringsList allowed_int_potential_names;  // Allowed interaction potential names
    StringsList allowed_ext_potential_names;  // Allowed external potential names
    StringsList allowed_coord_init_methods;   // Allowed coordinate initialization methods
    StringsList allowed_vel_init_methods;     // Allowed velocity initialization methods
    StringsList allowed_propagators;          // Allowed time propagation algorithms
    StringsList allowed_thermostats;          // Allowed thermostat names

    /// @todo Consider refactoring this method to a more general utility class (maybe to units)?
    static StringsList splitString(const std::string& line, char delimiter = ' ');
    static bool labelInArray(const std::string& label, const StringsList& arr);
};
