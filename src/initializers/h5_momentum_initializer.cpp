#ifdef USE_HDF5

#include "initializers/h5_momentum_initializer.h"
#include "initializers/h5_data_loader.h"
#include "core/system_state.h"

#include <format>

H5MomentumInitializer::H5MomentumInitializer(
    const std::string& filename,
    int first_idx,
    const std::string& init_vel_unit,
    long init_vel_frame,
    FrameSelectionMode init_vel_frame_mode,
    const std::shared_ptr<SystemState>& state,
    double mass
)
    : MomentumInitializer(state, mass),
    m_filename(filename),
    m_first_idx(first_idx),
    m_init_vel_unit(init_vel_unit),
    m_init_vel_frame(init_vel_frame),
    m_init_vel_frame_mode(init_vel_frame_mode) {
}

void H5MomentumInitializer::initialize() {
    const std::string h5_filename = std::vformat(
        m_filename, std::make_format_args(m_first_idx));

    H5DataLoader::loadFromFile(
        h5_filename,
        "velocities",
        m_init_vel_unit,
        "velocity",
        m_init_vel_frame,
        m_init_vel_frame_mode,
        m_state->momenta,
        m_mass
    );
}

#endif // USE_HDF5