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
}

class Params {
public:
	INIReader reader;

	using sys_types = std::variant<int, long, double, bool, std::string, std::vector<double>>;
	using sim_types = std::variant<int, long, double, bool, std::string, unsigned int>;
	using numbers_or_string = std::variant<int, long, double, std::string>;

	// Map holding the simulation parameters
	std::unordered_map<std::string, sim_types> sim;
	// Map holding the (physical) system parameters
	std::unordered_map<std::string, sys_types> sys;
	// Vector holding the interaction potential parameters
	std::vector<std::string> estimators;
	// Map holding the interaction potential parameters
	std::unordered_map<std::string, numbers_or_string> interaction_pot;
	// Map holding the interaction potential parameters
	std::unordered_map<std::string, numbers_or_string> external_pot;
	// Map holding the output parameters
	std::unordered_map<std::string, std::variant<int, double, bool>> out;

	Params(std::string filename);

	void setSimParam(std::string name, sim_types value);
	void setSysParam(std::string name, sys_types value);

	std::string getSysParam(std::string name);

	std::pair<double, std::string> parseQuantity(const std::string& input);
	double getQuantity(const std::string& family, const std::string& input);

	bool parseTokenParentheses(const std::string& input, std::string& token, std::string& value);

private:
	std::vector<std::string> allowed_estimators;           // Allowed estimators
	std::vector<std::string> interaction_potential_names;  // Allowed interaction potential names
	std::vector<std::string> external_potential_names;     // Allowed external potential names
	std::vector<std::string> allowed_coord_init_methods;   // Allowed coordinate initialization methods
	std::vector<std::string> allowed_vel_init_methods;     // Allowed velocity initialization methods
	std::vector<std::string> allowed_propagators;          // Allowed time propagation algorithms

	std::vector<std::string> splitString(const std::string& line, char delimiter = ' ');

	bool labelInArray(const std::string& label, const std::vector<std::string>& arr);
};