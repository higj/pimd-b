#include "initializers/xyz_position_initializer.h"
#include "initializers/xyz_data_loader.h"

XyzPositionInitializer::XyzPositionInitializer(
    const std::string& filename,
    int first_idx,
    const std::string& init_pos_unit,
    long init_pos_frame,
    XyzFrameSelectionMode init_pos_frame_mode,
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

void XyzPositionInitializer::initialize() {
    const std::string xyz_filename = std::vformat(
        m_filename,
        std::make_format_args(m_first_idx)
    );

    XyzDataLoader::loadFromFile(
        xyz_filename,
        m_init_pos_unit,
        "length",
        m_init_pos_frame,
        m_init_pos_frame_mode,
        *m_coord
    );
}