#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <format>

struct SimulationConfig;

/**
 * @brief Prints a summary of the simulation parameters at the end of the simulation.
 */
class SimulationReport {
public:
    explicit SimulationReport(const SimulationConfig& config, long steps);

    void writeReport(double wall_time);

private:
    // Pre-formatted report lines ready for output
    std::vector<std::string> m_parameter_lines;
    std::vector<std::string> m_feature_lines;

    // Data needed for runtime calculations
    long m_steps;
    bool m_should_skip_report;

    void initializeParameterLines(const SimulationConfig& config);
    void initializeFeatureLines();

    static std::string quantityInUserUnits(const double& value, const std::string& family, const std::string& key, const SimulationConfig& config);

    template <typename T>
    std::string formattedReportLine(const std::string& property_name, const T& value) {
        return std::format("{:<40}\t:\t{}\n", property_name, value);
    }
};