#include "initializers/xyz_momentum_initializer.h"
#include "initializers/xyz_data_loader.h"
#include "core/system_state.h"

#include <utility>

XyzMomentumInitializer::XyzMomentumInitializer(
    const std::string& filename,
    int first_idx,
    const std::string& init_vel_unit,
    long init_vel_frame,
    XyzFrameSelectionMode init_vel_frame_mode,
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

void XyzMomentumInitializer::initialize() {
    // Format filename with the current bead index (vformat returns literal if no fields)
    const std::string vel_filename = std::vformat(
        m_filename,
        std::make_format_args(m_first_idx)
    );

    XyzDataLoader::loadFromFile(
        vel_filename,
        m_init_vel_unit,
        "velocity",
        m_init_vel_frame,
        m_init_vel_frame_mode,
        m_state->momenta,
        m_mass
    );
}