#pragma once

#include "initializers/position_initializer.h"

#include <memory>

class RandomGenerators;

/**
 * @brief Samples random positions from a uniform distribution in the interval [-L/2, L/2],
 * for each particle along each axis.
 */
class RandomPositionInitializer final : public PositionInitializer {
public:
    explicit RandomPositionInitializer(const std::shared_ptr<RandomGenerators>& rng, const std::shared_ptr<dVec>& coords, double box_size);
    void initialize() override;
private:
    std::shared_ptr<RandomGenerators> m_rng;
};