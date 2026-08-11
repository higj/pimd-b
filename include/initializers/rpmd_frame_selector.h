#pragma once

#include <vector>
#include <string>
#include <memory>

struct SimulationConfig;

/**
 * @brief Handles RPMD frame selection from XYZ trajectory files.
 *
 * Computes the exact frame indices to use for each RPMD run based on:
 * - Total available frames in the NVT trajectory
 * - Number of independent RPMD runs to perform
 * - Equilibration discard fraction
 *
 * This class encapsulates XYZ file logic and should be used during Simulation
 * construction to prepare frame indices for RPMD runs.
 */
class RpmdFrameSelector {
public:
    /**
     * @brief Constructs the RPMD frame selector and computes frame indices.
     *
     * @param config Simulation configuration containing RPMD and XYZ file settings
     * @throws std::invalid_argument if RPMD is not enabled or config is invalid
     * @throws std::runtime_error if XYZ files cannot be read or contain no valid frames
     */
    explicit RpmdFrameSelector(const std::shared_ptr<SimulationConfig>& config);

    /**
     * @brief Returns the pre-computed frame indices for all RPMD runs.
     *
     * @return Vector of frame indices, one per RPMD run
     */
    [[nodiscard]] const std::vector<long>& getRpmdFrameIndices() const { return m_rpmd_frames; }

    /**
     * @brief Returns the total number of frames available in the XYZ file.
     *
     * @return Total frame count
     */
    [[nodiscard]] long getNumFrames() const { return m_num_frames; }

private:
    std::vector<long> m_rpmd_frames;
    long m_num_frames;

    /**
     * @brief Queries an XYZ file to determine the number of frames it contains.
     *
     * @param filename The XYZ file to query (may contain format placeholders for beads)
     * @param first_idx The bead index for filename formatting
     * @return Total number of frames in the file
     * @throws std::runtime_error if file cannot be opened or parsed
     */
    static long countFramesInXyzFile(const std::string& filename, int first_idx);

    /**
     * @brief Computes frame indices based on total frames, num_runs, and discard fraction.
     *
     * @param num_frames Total available frames
     * @param num_runs Number of RPMD runs
     * @param nvt_discard_frac Equilibration discard fraction
     * @return Vector of computed frame indices
     */
    static std::vector<long> computeFrameIndices(
        const long num_frames,
        const int num_runs,
        const double nvt_discard_frac
    );
};