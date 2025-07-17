#pragma once

#include "propagators/propagator.h"
#include "propagators/normal_modes/normal_modes.h"
#include "strategies/propagators/normal_modes_momenta_strategy.h"

struct BeadContext;
class BosonicExchangeBase;

class NormalModesPropagator final : public Propagator {
public:
    NormalModesPropagator(
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<ForceManager>& force_mgr,
        const BoxContext& box_ctx,
        const SpringContext& spring_ctx,
        double mass,
        double dt,
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
        const BeadContext& bead_ctx,
        bool bosonic,
        const std::shared_ptr<NormalModes>& normal_modes
    );
    ~NormalModesPropagator() override = default;

    void step() override;

private:
    std::shared_ptr<NormalModes> m_normal_modes;
    double m_freq, m_c, m_s, m_omega;
    std::unique_ptr<NormalModesMomentaStrategy> m_nm_momenta_strategy;

    void momentaExternalForces() const;
};
