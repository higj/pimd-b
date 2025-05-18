#pragma once

//#include "common.h"
#include "core/system_state.h"

class MomentumInitializer {
public:
    explicit MomentumInitializer(const std::shared_ptr<SystemState>& state, double mass)
        : m_state(state), m_natoms(m_state->getNumAtoms()), m_mass(mass) {}
    virtual ~MomentumInitializer() = default;
    virtual void initialize() = 0;
protected:
    std::shared_ptr<SystemState> m_state;
    int m_natoms;
    double m_mass;
};
