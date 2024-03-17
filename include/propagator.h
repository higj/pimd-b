#pragma once

class Simulation;

class Propagator {
public:
    explicit Propagator(Simulation& _sim);
    ~Propagator() = default;
    
    virtual void step() = 0;
protected:
    Simulation& sim; // Reference to the simulation object
};

/* -------------------------------- */

class VelocityVerletPropagator : public Propagator {
public:
    VelocityVerletPropagator(Simulation& _sim);
    
    void step() override;
};

/* -------------------------------- */

class NormalModesPropagator : public Propagator {
public:
    NormalModesPropagator(Simulation& _sim);
    
    void step() override;
};