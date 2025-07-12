#pragma once

#include "propagators/propagator.h"
#include "propagators/normal_modes/normal_modes.h"

class NormalModesPropagator final : public Propagator {
public:
    NormalModesPropagator(
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<ForceManager>& force_mgr,
        const std::shared_ptr<ExchangeState>& exchange_state,
        const SpringContext& spring_ctx,
        double mass,
        double dt,
        const std::shared_ptr<NormalModes>& normal_modes
    );
    ~NormalModesPropagator() override = default;

    void step() override;

private:
    std::shared_ptr<NormalModes> m_normal_modes;
    double m_freq, m_c, m_s, m_omega;

    void momentaExternalForces() const;
};
