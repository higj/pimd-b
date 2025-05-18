/*
 * This implementation is adapted from a similar code in i-PI.
 * Original source: https://github.com/i-pi/i-pi/blob/main/ipi/utils/units.py
 * License: MIT
 *
 * Copyright (c) [2014-2015] [i-PI developers]
 * This adaptation follows the terms of the original license.
 */

#pragma once

#include <unordered_map>
#include <string>

// Internally this simulation uses the atomic units.
// https://en.wikipedia.org/wiki/Hartree_atomic_units

namespace Constants {
    constexpr double kB = 1.0;         // Boltzmann constant (in Hartrees per atomic unit of temperature)
    constexpr double hbar = 1.0;       // Reduced Planck's constant
    constexpr double amu = 1822.8885;  // Atomic mass unit (in terms of electron masses)
}

/// @todo Add the option to detect short prefixes and units (e.g., "fs" or "femtosec" in addition to "femtoseconds")
namespace Units {
    // Metric prefixes
    const std::unordered_map<std::string, double> UnitPrefix = {
        {"", 1.0},
        {"yotta", 1e24},
        {"zetta", 1e21},
        {"exa", 1e18},
        {"peta", 1e15},
        {"tera", 1e12},
        {"giga", 1e9},
        {"mega", 1e6},
        {"kilo", 1e3},
        {"hecto", 1e2},
        {"deci", 1e-1},
        {"centi", 1e-2},
        {"milli", 1e-3},
        {"micro", 1e-6},
        {"nano", 1e-9},
        {"pico", 1e-12},
        {"femto", 1e-15},
        {"atto", 1e-18},
        {"zepto", 1e-21},
        {"yocto", 1e-24}
    };

    // Map of conversion factors from input unit to atomic units (which are used internally)
    const std::unordered_map<std::string, std::unordered_map<std::string, double>> UnitMap = {
        {"undefined", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00}
        }},
        {"energy", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
            {"electronvolt", 0.036749326},
            {"j/mol", 0.00000038087989},
            {"cal/mol", 0.0000015946679},
            {"kelvin", 3.1668152e-06},
        }},
        {"temperature", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
            {"kelvin", 3.1668152e-06},
        }},
        {"time", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
            {"second", 4.1341373e16},
        }},
        {"frequency", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
            {"inversecm", 4.5563353e-06},
            {"hertz*rad", 2.4188843e-17},
            {"hertz", 1.5198298e-16},
        }},
        {"ms-momentum", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
        }},
        {"length", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
            {"angstrom", 1.8897261},
            {"meter", 1.8897261e10},
            {"radian", 1.00},
            {"degree", 0.017453292519943295},
        }},
        {"volume", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
            {"angstrom3", 6.748334231},
        }},
        {"velocity", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
            {"angstrom/ps", 4.5710289e-5},
            {"m/s", 4.5710289e-7},  // WolframAlpha: 1/((Bohr radius)*(Hartree)/hbar)
        }},
        {"momentum", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
        }},
        {"mass", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
            {"dalton", 1.00 * Constants::amu},
            {"amu", 1.00 * Constants::amu},
            {"electronmass", 1.00},
        }},
        {"pressure", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
            {"bar", 3.398827377e-9},
            {"atmosphere", 3.44386184e-9},
            {"pascal", 3.398827377e-14},
            {"ev/ang3", 0.0054456877},
        }},
        {"density", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
            {"g/cm3", 162.67263},
        }},
        {"force", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
            {"newton", 12137805},
            {"ev/ang", 0.019446904},
        }},
        {"hessian", {
            {"", 1.00},
            {"automatic", 1.00},
            {"atomic_unit", 1.00},
            {"ev/ang^2", 0.010290858},
        }}
    };

    // Separates base unit from the metric prefix (if it exists)
    std::tuple<std::string, std::string> separatePrefixUnit(const std::string& unit);

    // Converts quantity from given units to internal (atomic) units
    double convertToInternal(const std::string& family, const std::string& unit, double number);

    // Converts the internal (atomic) units to the desired output units
    double convertToUser(const std::string& family, const std::string& unit, double number);

    // Parse a string containing a numerical value and a unit
    std::pair<double, std::string> parseQuantity(const std::string& input);

    /**
     * Accepts a string of the format "<number> <unit>" and returns the numerical
     * value of the quantity in internal (atomic) units.
     *
     * @param family The family of units to which the quantity belongs.
     * @param input The string containing the numerical value and the unit.
     * @return The numerical value of the quantity in internal units.
     */
    double getQuantity(const std::string& family, const std::string& input);
}
