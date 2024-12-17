#pragma once

#include <vector>
#include "mpi.h"
#include "common.h"

#ifndef CALC_ETA
#define CALC_ETA false
#endif

class Simulation;

class Thermostat {
public:
    explicit Thermostat(Simulation& _sim);
    virtual ~Thermostat() = default;

    virtual void step();
protected:
    Simulation& sim; // Reference to the simulation object
};

/* -------------------------------- */

class LangevinThermostat : public Thermostat {
public:
    LangevinThermostat(Simulation& _sim);
    ~LangevinThermostat() override = default;

    void step() override;
private:
    double a, b;
};

/* -------------------------------- */

class NoseHooverThermostat : public Thermostat {
public:
    NoseHooverThermostat(Simulation& _sim, int _nchains);
    ~NoseHooverThermostat() override = default;

    void step() override;
protected:
    double singleChainStep(const double& current_energy, const int& index);
    int nchains;
    double Q1, Qi;
    double dt2, dt4, dt8;
    double required_energy;
    std::vector<double> eta_dot, eta_dot_dot; 
#if CALC_ETA
    std::vector<double> eta; 
#endif
};

/* -------------------------------- */

class NoseHooverNpThermostat : public NoseHooverThermostat {
public:
    NoseHooverNpThermostat(Simulation& _sim, int _nchains);
    ~NoseHooverNpThermostat() override = default;

    void step() override;
};
/* -------------------------------- */

class NoseHooverNpDimThermostat : public NoseHooverThermostat {
public:
    NoseHooverNpDimThermostat(Simulation& _sim, int _nchains);
    ~NoseHooverNpDimThermostat() override = default;

    void step() override;
};
/* -------------------------------- */
