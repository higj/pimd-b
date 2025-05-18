#pragma once

#include "contexts/propagator_context.h"

class Propagator {
public:
    explicit Propagator(const PropagatorContext& context);
    virtual ~Propagator() = default;
    
    virtual void step() = 0;
    void momentStep() const;
    void coordsStep() const;

protected:
    PropagatorContext m_context;
};