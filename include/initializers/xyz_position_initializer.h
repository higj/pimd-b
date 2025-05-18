#pragma once

#include "initializers/position_initializer.h"

class XyzPositionInitializer final : public PositionInitializer {
public:
    XyzPositionInitializer(const std::string& filename, int first_idx, const std::shared_ptr<dVec>& coord, double box_size);

    /**
     * Initializes positions from user-provided files.
     */
    void initialize() override;

private:
    std::string m_filename;
    int m_first_idx;

    /**
     * Load positions from an .xyz file.
     *
     * @param pos_filename Filename of the positions file.
     * @param destination Destination vector to store the loaded momenta.
     */
    void loadFromFile(const std::string& pos_filename, dVec& destination) const;
};