#pragma once

#include "distinguishable_normal_modes_momenta_strategy.h"

class BosonicExchangeBase;

class BosonicNormalModesMomentaStrategy : public DistinguishableNormalModesMomentaStrategy {
public:
    BosonicNormalModesMomentaStrategy(
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange
    );

    void momentaExternalForces(
        const std::shared_ptr<SystemState>& state,
        const SpringContext& spring_ctx,
        double dt
    ) override;

private:
    std::shared_ptr<BosonicExchangeBase> m_bosonic_exchange;
};
