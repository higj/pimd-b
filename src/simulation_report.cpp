#include "simulation_report.h"
#include "core/simulation_config.h"
#include "output_paths.h"
#include "units.h"

#include <chrono>

SimulationReport::SimulationReport(const SimulationConfig& config, long steps)
    : m_steps(steps), m_should_skip_report(config.this_bead != 0) {
    if (!m_should_skip_report) {
        initializeParameterLines(config);
        initializeFeatureLines();
    }
}

void SimulationReport::writeReport(double wall_time) {
    if (m_should_skip_report)
        return;

    std::ofstream report_file;
    report_file.open(std::format("{}/report.txt", Output::FOLDER_NAME), std::ios::out | std::ios::app);

    // Write parameters section
    report_file << "---------\nParameters\n---------\n";
    for (const auto& line : m_parameter_lines) {
        report_file << line;
    }

    // Write features section
    report_file << "---------\nFeatures\n---------\n";
    for (const auto& line : m_feature_lines) {
        report_file << line;
    }

    // Write runtime information (calculated at write time)
    report_file << "---------\n";
    report_file << formattedReportLine("Wall time", std::format("{:%T}",
        std::chrono::duration<double>(wall_time)));
    report_file << formattedReportLine("Wall time per step (sec)",
        std::format("{:.5e}", wall_time / m_steps));

    report_file.close();
}

std::string SimulationReport::quantityInUserUnits(
    const double& value,
    const std::string& family, 
    const std::string& key, 
    const SimulationConfig& config)
{
    const std::string out_unit = config.units_list.at(key);
    double out_value = Units::convertToUser(family, out_unit, value);
    return std::format("{} {}", out_value, out_unit);
}

void SimulationReport::initializeParameterLines(const SimulationConfig& config) {
    // Statistics and algorithm information
    if (config.bosonic) {
        m_parameter_lines.push_back(formattedReportLine("Statistics", "Bosonic"));

        std::string bosonic_alg_name = "Feldman-Hirshberg";
#if FACTORIAL_BOSONIC_ALGORITHM
        bosonic_alg_name = "Naive";
#endif
        m_parameter_lines.push_back(formattedReportLine("Bosonic algorithm", bosonic_alg_name));
    } else {
        m_parameter_lines.push_back(formattedReportLine("Statistics", "Boltzmannonic"));
    }

    /// TODO: Add thermostat

    // Configuration parameters
    m_parameter_lines.push_back(formattedReportLine("Time propagation algorithm", config.propagator_type));
    m_parameter_lines.push_back(
        formattedReportLine(
            "Time step", 
            quantityInUserUnits(config.dt, "time", "dt", config)
        )
    );
    m_parameter_lines.push_back(formattedReportLine("Periodic boundary conditions", config.pbc));
    m_parameter_lines.push_back(formattedReportLine("Dimension", NDIM));
    m_parameter_lines.push_back(formattedReportLine("Seed", config.seed));
    m_parameter_lines.push_back(formattedReportLine("Coordinate initialization method", config.init_pos_type));
    m_parameter_lines.push_back(formattedReportLine("Momentum initialization method", config.init_vel_type));
    m_parameter_lines.push_back(formattedReportLine("Number of atoms", config.natoms));
    m_parameter_lines.push_back(formattedReportLine("Number of beads", config.nbeads));

    // Unit conversions
    /*
    double out_temperature = Units::convertToUser("temperature", "kelvin", config.temperature);
    m_parameter_lines.push_back(formattedReportLine("Temperature", std::format("{} kelvin", out_temperature)));

    double out_sys_size = Units::convertToUser("length", "angstrom", config.box_size);
    m_parameter_lines.push_back(formattedReportLine("Linear size of the system", std::format("{} angstroms", out_sys_size)));

    double out_mass = Units::convertToUser("mass", "dalton", config.mass);
    m_parameter_lines.push_back(formattedReportLine("Mass", std::format("{} amu", out_mass)));
    */
    m_parameter_lines.push_back(
        formattedReportLine(
            "Temperature",
            quantityInUserUnits(config.temperature, "temperature", "temperature", config)
        )
    );

    m_parameter_lines.push_back(
        formattedReportLine(
            "Linear size of the system",
            quantityInUserUnits(config.box_size, "length", "size", config)
        )
    );

    m_parameter_lines.push_back(
        formattedReportLine(
            "Mass",
            quantityInUserUnits(config.mass, "mass", "mass", config)
        )
    );

    m_parameter_lines.push_back(formattedReportLine("Total number of MD steps", config.steps));
    m_parameter_lines.push_back(formattedReportLine("Interaction potential name", config.int_potential_cfg.name()));
    m_parameter_lines.push_back(formattedReportLine("External potential name", config.ext_potential_cfg.name()));

    if (config.rpmd_config.enabled)
    {
        //m_parameter_lines.push_back(formattedReportLine("RPMD mode", "Enabled"));
        m_parameter_lines.push_back(formattedReportLine("RPMD number of runs", config.rpmd_config.num_runs));
        m_parameter_lines.push_back(formattedReportLine("RPMD discard fraction", config.rpmd_config.nvt_discard_frac));
    }
}

void SimulationReport::initializeFeatureLines() {
    m_feature_lines.push_back(formattedReportLine("Minimum image convention", MINIM));
    m_feature_lines.push_back(formattedReportLine("Wrapping of coordinates", WRAP));
    m_feature_lines.push_back(formattedReportLine("Tau convention", TAU_CONVENTION));
}
