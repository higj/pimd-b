#pragma once

#include <vector>

#include "thermostats/thermostat.h"

class Simulation;
class Coupling;

class NoseHooverThermostat : public Thermostat
{
public:
    NoseHooverThermostat(
        const ThermalContext& thermal_ctx,
        const NormalModesContext& nm_ctx,
        const std::shared_ptr<SystemState>& state,
        int nchains,
        double dt,
        double mass
    );
    ~NoseHooverThermostat() override = default;

    void momentaUpdate() override;

    /**
     * Returns the contribution of NHC to the conserved quantity.
     * Corresponds to the difference between the conserved quantity
     * and the classical Hamiltonian of the ring polymers.
     * The equations of motion conserve H + additionToH, where H
     * is the Hamiltonian of the physical system
     */
    double getAdditionToH() override;
protected:
    /**
     * Calculates the scaling factor of the momenta based on
     * Tuckerman et al. (2006) J. Phys. A: Math. Gen. 39 5629
     * and its implementation in LAMMPS.
     *
     * @param current_energy The current energy
     * @param index          The index of the first component in the chain
     * @return               The scaling factor
     */
    double singleChainStep(const double& current_energy, const int& index);

    int m_nbeads;
    int m_natoms;
    int m_nchains; // Number of components in each Nose-Hoover chain
    double m_dt; // Time step (TODO: dt needed only in the constructor; we later use only dt2, dt4 and dt8)
    double m_mass; // Mass of the particles

    double Q1, Qi; // The masses of the eta
    double dt2, dt4, dt8;
    double required_energy;
    // A measure of the required preserved energy (kT times the number of degrees of freedom associated with the chain)
    std::vector<double> eta_dot, eta_dot_dot; // The first and second derivatives of eta
    std::vector<double> eta;

    /**
     * Calculate the contribution of NHC to the conserved quantity.
     *
     * @param ndof  Number of degrees of freedom
     * @param index Index of the first component in the chain
     * @return      Energy contribution to the conserved quantity
     */
    [[nodiscard]] double singleChainGetAdditionToH(const int& ndof, const int& index) const;
};

/* -------------------------------- */

class NoseHooverNpThermostat final : public NoseHooverThermostat
{
public:
    NoseHooverNpThermostat(
        const ThermalContext& thermal_ctx,
        const NormalModesContext& nm_ctx,
        const std::shared_ptr<SystemState>& state,
        int nchains,
        double dt,
        double mass
    );
    ~NoseHooverNpThermostat() override = default;

    void momentaUpdate() override;
    double getAdditionToH() override;
    // The equations of motion conserve H + additionToH, where H is the Hamiltonian of the physical system
};

/* -------------------------------- */

class NoseHooverNpDimThermostat final : public NoseHooverThermostat
{
public:
    NoseHooverNpDimThermostat(
        const ThermalContext& thermal_ctx,
        const NormalModesContext& nm_ctx,
        const std::shared_ptr<SystemState>& state,
        int nchains,
        double dt,
        double mass
    );
    ~NoseHooverNpDimThermostat() override = default;

    void momentaUpdate() override;
    double getAdditionToH() override;
    // The equations of motion conserve H + additionToH, where H is the Hamiltonian of the physical system
};
