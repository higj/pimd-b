#include "initializers/rpmd_frame_selector.h"
#include "core/simulation_config.h"

#include <fstream>
#include <sstream>
#include <format>
#include <cmath>
#include <stdexcept>

RpmdFrameSelector::RpmdFrameSelector(const std::shared_ptr<SimulationConfig>& config) {
    if (!config) {
        throw std::invalid_argument("SimulationConfig cannot be null");
    }

    if (!config->rpmd_config.enabled) {
        throw std::invalid_argument("RPMD frame selector requires RPMD to be enabled");
    }

    if (config->init_pos_type != "xyz" || config->init_vel_type != "xyz") {
        throw std::invalid_argument("RPMD frame selector requires both position and velocity initialization files to be 'xyz'");
    }

    // Determine which file to query (position or velocity, both should have same frame count)
    const std::string& filename = config->init_pos_filename;
    const int first_idx = config->this_bead + config->init_pos_index_offset;

    // Query XYZ file for total frame count
    m_num_frames = countFramesInXyzFile(filename, first_idx);

    // Compute frame indices for all RPMD runs
    m_rpmd_frames = computeFrameIndices(
        m_num_frames,
        config->rpmd_config.num_runs,
        config->rpmd_config.nvt_discard_frac
    );
}

long RpmdFrameSelector::countFramesInXyzFile(const std::string& filename, int first_idx) {
    // Format the filename if it contains a placeholder for the bead index
    const std::string actual_filename = std::vformat(filename, std::make_format_args(first_idx));

    std::ifstream file(actual_filename);
    if (!file.is_open()) {
        throw std::runtime_error(
            std::format("Cannot open XYZ file for RPMD frame counting: {}", actual_filename)
        );
    }

    long frame_count = 0;
    std::string line;

    // XYZ format: first line of each frame contains the number of atoms
    // We count frames by counting valid lines that contain atom counts
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty() || line.find_first_not_of(" \t\n\r") == std::string::npos) {
            continue;
        }

        // Try to parse as atom count
        try {
            std::istringstream iss(line);
            int atom_count;
            if (iss >> atom_count && atom_count > 0) {
                // This looks like a valid frame header
                frame_count++;

                // Skip the next atom_count + 1 lines (comment line + atom data)
                for (int i = 0; i < atom_count + 1; ++i) {
                    if (!std::getline(file, line)) {
                        throw std::runtime_error(
                            std::format("Unexpected EOF in XYZ file: {}", actual_filename)
                        );
                    }
                }
            }
        } catch (const std::exception&) {
            // Continue on parse errors; XYZ format is forgiving
            continue;
        }
    }

    if (frame_count == 0) {
        throw std::runtime_error(
            std::format("No valid frames found in XYZ file: {}", actual_filename)
        );
    }

    return frame_count;
}

std::vector<long> RpmdFrameSelector::computeFrameIndices(
    const long num_frames,
    const int num_runs,
    const float nvt_discard_frac
)
{
    const long equilibration_frame = static_cast<long>(
        std::floor(nvt_discard_frac * num_frames)
        );

    std::vector<long> rpmd_frames;

    if (num_runs == 1) {
        rpmd_frames.push_back(equilibration_frame);
    } else {
        const long last_frame = num_frames - 1;

        for (int run = 0; run < num_runs; ++run) {
            const double fraction = static_cast<double>(run) / (num_runs - 1);
            const long frame = equilibration_frame + static_cast<long>(
                std::floor(fraction * (last_frame - equilibration_frame))
                );

            rpmd_frames.push_back(frame);
        }
    }

    return rpmd_frames;
}