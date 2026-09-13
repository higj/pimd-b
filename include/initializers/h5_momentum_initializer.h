#pragma once
#ifdef USE_HDF5

#include "core/simulation_config.h"
#include "initializers/momentum_initializer.h"

#include <memory>
#include <string>

/**
 * @brief Reads initial velocities from one frame of an HDF5 trajectory file
 *        and converts them to momenta.
 *
 * The HDF5 file must contain a 3D dataset "velocities" with shape
 * [n_frames, n_atoms, n_coords]. Each value is multiplied by mass after
 * unit conversion to produce momenta, mirroring XyzMomentumInitializer.
 */
class H5MomentumInitializer final : public MomentumInitializer {
public:
    H5MomentumInitializer(
        const std::string& filename,
        int first_idx,
        const std::string& init_vel_unit,
        long init_vel_frame,
        FrameSelectionMode init_vel_frame_mode,
        const std::shared_ptr<SystemState>& state,
        double mass
    );

    void initialize() override;

private:
    std::string m_filename;
    int m_first_idx;
    std::string m_init_vel_unit;
    long m_init_vel_frame;
    FrameSelectionMode m_init_vel_frame_mode;
};

#endif // USE_HDF5