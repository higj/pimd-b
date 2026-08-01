#pragma once

#include "core/simulation_config.h"
#include "initializers/momentum_initializer.h"

#include <memory>
#include <string>

class XyzMomentumInitializer final : public MomentumInitializer {
public:
    /**
     * @brief Creates an initializer that reads one frame from an XYZ trajectory for velocities.
     *
     * The filename may have one format replacement field, which is replaced
     * with first_idx so individual MPI beads can use separate input files.
     * init_vel_frame is either a zero-based frame index or a value from a
     * "Step <number>" comment, according to init_vel_frame_mode.
     *
     * @param filename XYZ file path (may include format field like "{}")
     * @param first_idx Index to substitute in filename format
     * @param init_vel_unit Unit of velocities in the XYZ file
     * @param init_vel_frame Frame index or step number to load
     * @param init_vel_frame_mode How to interpret init_vel_frame (Index or Step)
     * @param state Simulation system state
     * @param mass Particle mass (used for velocity to momentum conversion)
     */
    XyzMomentumInitializer(
        const std::string& filename,
        int first_idx,
        const std::string& init_vel_unit,
        long init_vel_frame,
        XyzFrameSelectionMode init_vel_frame_mode,
        const std::shared_ptr<SystemState>& state,
        double mass
    );

    /**
     * Initializes momenta from user-provided XYZ files.
     */
    void initialize() override;

private:
    std::string m_filename;
    int m_first_idx;
    std::string m_init_vel_unit;
    long m_init_vel_frame;
    XyzFrameSelectionMode m_init_vel_frame_mode;
};