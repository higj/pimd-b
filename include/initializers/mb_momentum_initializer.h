#pragma once

#include "initializers/momentum_initializer.h"

#include <memory>

class RandomGenerators;

class MaxwellBoltzmannMomentumInitializer final : public MomentumInitializer {
public:
    MaxwellBoltzmannMomentumInitializer(const std::shared_ptr<RandomGenerators>& rng, const std::shared_ptr<SystemState>& state, double mass, double thermo_beta);
    void initialize() override;

private:
    std::shared_ptr<RandomGenerators> m_rng;
    double m_thermo_beta;
    double sampleMaxwellBoltzmann(double thermo_beta, double mass) const;
};