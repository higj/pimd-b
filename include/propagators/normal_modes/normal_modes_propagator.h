#pragma once

#include "propagators/propagator.h"
#include "propagators/normal_modes/normal_modes.h"

class NormalModesPropagator final : public Propagator {
public:
    NormalModesPropagator(const PropagatorContext& context, const std::shared_ptr<NormalModes>& normal_modes);
    ~NormalModesPropagator() override = default;

    void step() override;

    std::shared_ptr<NormalModes> m_normal_modes;
private:
    double m_freq, m_c, m_s, m_omega;
    void momentaExternalForces() const;
};
