#pragma once

#include "initializers/momentum_initializer.h"

class ManualMomentumInitializer final : public MomentumInitializer {
public:
    ManualMomentumInitializer(const std::string& filename, int first_idx, const std::shared_ptr<SystemState>& state, double mass);

    /**
     * Initializes momenta from user-provided files.
     */
    void initialize() override;

private:
    std::string m_filename;
    int m_first_idx;

    /**
     * Load momenta from a file.
     *
     * @param vel_filename Filename of the velocity file.
     * @param destination Destination vector to store the loaded momenta.
     */
    void loadFromFile(const std::string& vel_filename, dVec& destination) const;
};