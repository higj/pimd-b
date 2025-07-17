#pragma once

#include "initializers/position_initializer.h"

#include <memory>

/**
 * @brief Arranges particles in a grid-like structure within the simulation box.
 */
class GridPositionInitializer final : public PositionInitializer {
public:
    explicit GridPositionInitializer(const std::shared_ptr<dVec>& coords, const BoxContext& box_ctx);
    void initialize() override;
};