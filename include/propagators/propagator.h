#pragma once

#include "contexts/spring_context.h"
#include <memory>

class SystemState;
class ForceManager;
class ExchangeState;

class Propagator
{
public:
    explicit Propagator(
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<ForceManager>& force_mgr,
        const std::shared_ptr<ExchangeState>& exchange_state,
        const SpringContext& spring_ctx,
        double mass,
        double dt
    );
    virtual ~Propagator() = default;

    virtual void step() = 0;

protected:
    // Propagators need access to the system state (including updateNeighboringCoordinates)
    std::shared_ptr<SystemState> m_state;
    // Propagators need access to the force field manager (including updateSpringForces)
    std::shared_ptr<ForceManager> m_force_mgr;
    // Propagators need access to bosonic exchange methods
    std::shared_ptr<ExchangeState> m_exchange_state;
    // Propagators need access to the spring context
    SpringContext m_spring_ctx;

    int m_natoms;
    int m_nbeads;
    int m_this_bead;

    double m_mass;
    double m_dt;
};