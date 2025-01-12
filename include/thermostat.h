#pragma once

#include <vector>
#include "common.h"
#include <memory>

class Simulation;
class Coupling;

class Thermostat {
public:
    explicit Thermostat(Simulation& _sim, bool normal_modes);
    virtual ~Thermostat() = default;
    void step();
    virtual void momentaUpdate();
    virtual double getAdditionToH();
protected:
    Simulation& sim;   // Reference to the simulation object
    std::unique_ptr<Coupling> coupling;
};

/* -------------------------------- */

class LangevinThermostat : public Thermostat {
public:
    LangevinThermostat(Simulation& _sim, bool normal_modes);
    ~LangevinThermostat() override = default;

    void momentaUpdate() override;

protected:
    double friction_coefficient, noise_coefficient;
};

/* -------------------------------- */

class NoseHooverThermostat : public Thermostat {
public:
    NoseHooverThermostat(Simulation& _sim, bool normal_modes, int _nchains);
    ~NoseHooverThermostat() override = default;

    void momentaUpdate() override;
    double getAdditionToH() override; // The equations of motion conserve H + additionToH, where H is the Hamiltonian of the physical system
protected:
    double singleChainStep(const double& current_energy, const int& index);
    int nchains; // Number of components in each Nose-Hoover chain
    double Q1, Qi; // The masses of the eta
    double dt2, dt4, dt8;
    double required_energy; // A measure of the required preserved energy (kT times the number of degrees of freedom associated with the chain)
    std::vector<double> eta_dot, eta_dot_dot; // The first and second derivatives of eta
    std::vector<double> eta; 
    double singleChainGetAdditionToH(const int& expected_energy, const int& index);
};

/* -------------------------------- */

class NoseHooverNpThermostat : public NoseHooverThermostat {
public:
    NoseHooverNpThermostat(Simulation& _sim, bool normal_modes, int _nchains);
    ~NoseHooverNpThermostat() override = default;

    void momentaUpdate() override;
    double getAdditionToH() override; // The equations of motion conserve H + additionToH, where H is the Hamiltonian of the physical system
};

/* -------------------------------- */

class NoseHooverNpDimThermostat : public NoseHooverThermostat {
public:
    NoseHooverNpDimThermostat(Simulation& _sim, bool normal_modes, int _nchains);
    ~NoseHooverNpDimThermostat() override = default;

    void momentaUpdate() override;
    double getAdditionToH() override; // The equations of motion conserve H + additionToH, where H is the Hamiltonian of the physical system
};
