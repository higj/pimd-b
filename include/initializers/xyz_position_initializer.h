#pragma once

#include "core/simulation_config.h"
#include "initializers/position_initializer.h"

#include <memory>
#include <string>

class XyzPositionInitializer final : public PositionInitializer {
public:
    /**
     * @brief Creates an initializer that reads one frame from an XYZ trajectory.
     *
     * The filename may have one format replacement field, which is replaced
     * with first_idx so individual MPI beads can use separate input files.
     * init_pos_frame is either a zero-based frame index or a value from a
     * "Step <number>" comment, according to init_pos_frame_mode.
     */
    XyzPositionInitializer(
        const std::string& filename,
        int first_idx,
        const std::string& init_pos_unit,
        long init_pos_frame,
        XyzFrameSelectionMode init_pos_frame_mode,
        const std::shared_ptr<VecArray>& coord,
        const BoxContext& box_ctx
    );

    /**
     * Initializes positions from user-provided files.
     */
    void initialize() override;

private:
    std::string m_filename;
    int m_first_idx;
    std::string m_init_pos_unit;
    long m_init_pos_frame;
    XyzFrameSelectionMode m_init_pos_frame_mode;
};