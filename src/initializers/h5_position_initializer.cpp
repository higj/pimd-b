#ifdef USE_HDF5

#include "initializers/h5_position_initializer.h"
#include "initializers/h5_data_loader.h"

#include <format>

H5PositionInitializer::H5PositionInitializer(
    const std::string& filename,
    int first_idx,
    const std::string& init_pos_unit,
    long init_pos_frame,
    FrameSelectionMode init_pos_frame_mode,
    const std::shared_ptr<VecArray>& coord,
    const BoxContext& box_ctx
)
    : PositionInitializer(coord, box_ctx),
    m_filename(filename),
    m_first_idx(first_idx),
    m_init_pos_unit(init_pos_unit),
    m_init_pos_frame(init_pos_frame),
    m_init_pos_frame_mode(init_pos_frame_mode) {
}

void H5PositionInitializer::initialize() {
    const std::string h5_filename = std::vformat(
        m_filename, std::make_format_args(m_first_idx));

    H5DataLoader::loadFromFile(
        h5_filename,
        "positions",
        m_init_pos_unit,
        "length",
        m_init_pos_frame,
        m_init_pos_frame_mode,
        *m_coord
    );
}

#endif // USE_HDF5