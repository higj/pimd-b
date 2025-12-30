#pragma once

#include "initializers/position_initializer.h"

#include <memory>

class XyzPositionInitializer final : public PositionInitializer {
public:
    XyzPositionInitializer(const std::string& filename, int first_idx, const std::string& init_pos_unit, const std::shared_ptr<dVec>& coord, const BoxContext& box_ctx);

    /**
     * Initializes positions from user-provided files.
     */
    void initialize() override;

private:
    std::string m_filename;
    int m_first_idx;
    std::string m_init_pos_unit;

    /**
     * Load positions from an .xyz file.
     *
     * @param xyz_filename Filename of the positions file.
     * @param init_pos_unit Unit of the positions in the file.
     * @param destination Destination vector to store the loaded momenta.
     */
    static void loadFromFile(const std::string& xyz_filename, const std::string& init_pos_unit, dVec& destination);
};