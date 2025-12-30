#pragma once

#include "initializers/momentum_initializer.h"

#include <memory>

class ManualMomentumInitializer final : public MomentumInitializer {
public:
    ManualMomentumInitializer(const std::string& filename, int first_idx, const std::string& init_vel_unit, const std::shared_ptr<SystemState>& state, double mass);

    /**
     * Initializes momenta from user-provided files.
     */
    void initialize() override;

private:
    std::string m_filename;
    int m_first_idx;
    std::string m_init_vel_unit;

    /**
     * Load momenta from a file.
     *
     * @param vel_filename Filename of the velocity file.
     * @param init_vel_unit Unit of the velocities in the file.
     * @param destination Destination vector to store the loaded momenta.
     */
    void loadFromFile(const std::string& vel_filename, const std::string& init_vel_unit, dVec& destination) const;
};