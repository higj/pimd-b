#pragma once
#ifdef USE_HDF5

#include "core/simulation_config.h"
#include "initializers/position_initializer.h"

#include <memory>
#include <string>

/**
 * @brief Reads initial positions from one frame of an HDF5 trajectory file.
 *
 * The filename may contain one format field replaced by first_idx (same
 * per-bead convention as XyzPositionInitializer). Frame selection is by
 * ordinal index or by matching step number, controlled by init_pos_frame_mode.
 *
 * The HDF5 file must contain a 3D dataset "positions" with shape
 * [n_frames, n_atoms, n_coords].
 */
class H5PositionInitializer final : public PositionInitializer {
public:
    H5PositionInitializer(
        const std::string& filename,
        int first_idx,
        const std::string& init_pos_unit,
        long init_pos_frame,
        FrameSelectionMode init_pos_frame_mode,
        const std::shared_ptr<VecArray>& coord,
        const BoxContext& box_ctx
    );

    void initialize() override;

private:
    std::string m_filename;
    int m_first_idx;
    std::string m_init_pos_unit;
    long m_init_pos_frame;
    FrameSelectionMode m_init_pos_frame_mode;
};

#endif // USE_HDF5