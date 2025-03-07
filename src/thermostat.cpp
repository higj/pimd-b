#include "thermostat.h"
#include "simulation.h"
#include "common.h"
#include "normal_modes.h"
#include "thermostat_coupling.h"

#include <numbers>

Thermostat::Thermostat(Simulation& _sim, bool normal_modes) : sim(_sim) {
    // Choose coupling (Cartesian coords or normal modes of distinguishable ring polymers)
    if (normal_modes) {
        coupling = std::make_unique<NMCoupling>(_sim);
    } else {
        coupling = std::make_unique<CartesianCoupling>(_sim);
    }
}

// This is the step function of a general thermostat, called in the simulation's run loop
void Thermostat::step() {
    coupling->mpiCommunication();
    momentaUpdate();
    coupling->updateCoupledMomenta();
}

// This is an update of the momenta within the thermostat step, unique for each thermostat
void Thermostat::momentaUpdate() {}

double Thermostat::getAdditionToH() {
    return 0;
}

/* -------------------------------- */

LangevinThermostat::LangevinThermostat(Simulation& _sim, bool normal_modes) : Thermostat(_sim, normal_modes) {
    friction_coefficient = exp(-0.5 * sim.gamma * sim.dt);
    noise_coefficient = sqrt((1 - friction_coefficient * friction_coefficient) * sim.mass / sim.thermo_beta);
}

void LangevinThermostat::momentaUpdate() {
    //std::normal_distribution<double> normal; // At the moment we use the Marsaglia generator

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double noise = sim.mars_gen->gaussian(); // LAMMPS pimd/langevin random numbers
            double momentum_for_calc = coupling->getMomentumForCalc(ptcl_idx, axis);
            double& momentum_for_update = coupling->getMomentumForUpdate(ptcl_idx, axis);
            // Perturb the momenta with a Langevin thermostat
            momentum_for_update = friction_coefficient * momentum_for_calc + noise_coefficient * noise;
        }
    }
}

/* -------------------------------- */

NoseHooverThermostat::NoseHooverThermostat(Simulation& _sim, bool normal_modes, int _nchains) : Thermostat(
    _sim, normal_modes) {
    nchains = _nchains;

    // Q1 and Qi (where i>1) are the masses of the extended variables. The choice
    // is based on Sec. 2.5 in Martyna, G. J. et al. (1996), Mol. Phys., 87(5), pp. 1117-1157.
    // For the fluctuation frequency, we choose the spring frequency of the ring polymers.
#if IPI_CONVENTION
    /// TODO: Explain why there is no division by nbeads
    Qi = Constants::hbar * Constants::hbar * sim.beta;
#else
    Qi = Constants::hbar * Constants::hbar * sim.beta / sim.nbeads;
#endif
    // Q1/Qi should be equal to the number of degrees of freedom (NDOF).
    // In the absence of constraints, NDOF=NDIM*natoms.
    Q1 = NDIM * sim.natoms * Qi;

    eta = std::vector<double>(nchains, 0.0);
    eta_dot = std::vector<double>(nchains, 0.0);
    eta_dot_dot = std::vector<double>(nchains, 0.0);

    dt2 = 0.5 * sim.dt;
    dt4 = 0.25 * sim.dt;
    dt8 = 0.125 * sim.dt;

    required_energy = NDIM * sim.natoms / sim.thermo_beta;
}

// CR: Add a comment explaining what the function does
double NoseHooverThermostat::getAdditionToH() {
    return singleChainGetAdditionToH(NDIM * sim.natoms, 0);
}

/**
 * Calculate the contribution of NHC to the conserved quantity.
 * Corresponds to the difference between the conserved quantity
 * and the classical Hamiltonian of the ring polymers.
 *
 * @param ndof Number of degrees of freedom
 * @param index Index of the first component in the chain.
 * @return Energy contribution to the conserved quantity.
 */
double NoseHooverThermostat::singleChainGetAdditionToH(const int& ndof, const int& index) {
    // First, we add the kinetic energy associated with the thermostat.
    // Note that the momentum of the thermostat is given by Qi*eta_dot.
    double addition_to_H = 0.5 * Q1 * eta_dot[index] * eta_dot[index];

    // Then, we add the other term that is proportional to kT*eta
    addition_to_H += ndof * eta[index] / sim.thermo_beta;

    // Similar calculation, but for the rest of the thermostats
    for (int i = 1; i < nchains; i++) {
        addition_to_H += 0.5 * Qi * eta_dot[index + i] * eta_dot[index + i];
        addition_to_H += eta[index + i] / sim.thermo_beta;
    }

    return addition_to_H;
}

void NoseHooverThermostat::momentaUpdate() {
    // Calculates the current energy in the system
    double current_energy = 0.0;
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            const double momentum_for_calc = coupling->getMomentumForCalc(ptcl_idx, axis);
            current_energy += momentum_for_calc * momentum_for_calc;
        }
    }
    current_energy *= 1 / sim.mass;

    // Obtain the scaling factor
    double scale = singleChainStep(current_energy, 0);

    // Rescale the momenta
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double momentum_for_calc = coupling->getMomentumForCalc(ptcl_idx, axis);
            double& momentum_for_update = coupling->getMomentumForUpdate(ptcl_idx, axis);
            momentum_for_update = momentum_for_calc * scale;
        }
    }
}

/**
 * Calculates the scaling factor of the momenta based on
 * Tuckerman et al. (2006) J. Phys. A: Math. Gen. 39 5629
 * and its implementation in LAMMPS.
 * 
 * @param current_energy The current energy
 * @param index          The index of the first component in the chain
 * @return               The scaling factor 
 */
double NoseHooverThermostat::singleChainStep(const double& current_energy, const int& index) {
    double exp_factor = 0.0;

    // Update the second derivative of eta for the first component
    eta_dot_dot[index] = (current_energy - required_energy) / Q1;

    // Update the first derivative of eta for the last component
    eta_dot[nchains - 1 + index] += eta_dot_dot[nchains - 1 + index] * dt4;

    // Update the first derivative of eta for all others components
    for (int i = nchains - 2; i >= 0; i--) {
        exp_factor = exp(-dt8 * eta_dot[i + 1 + index]);
        eta_dot[i + index] *= exp_factor;
        eta_dot[i + index] += eta_dot_dot[i + index] * dt4;
        eta_dot[i + index] *= exp_factor;
    }

    // Calculate the rescaling factor
    double scale = exp(-dt2 * eta_dot[index]);
    for (int i = 0; i < nchains; i++)
        eta[i + index] += dt2 * eta_dot[i + index];

    // Update the derivatives of eta for the first component
    eta_dot_dot[index] = (current_energy * scale * scale - required_energy) / Q1;
    eta_dot[index] *= exp_factor;
    eta_dot[index] += eta_dot_dot[index] * dt4;
    eta_dot[index] *= exp_factor;

    // Update the derivatives of eta for all other components except the last
    double Q_former = Q1;
    for (int i = 1; i < nchains - 1; i++) {
        exp_factor = exp(-dt8 * eta_dot[i + 1 + index]);
        eta_dot[i + index] *= exp_factor;
        eta_dot_dot[i + index] = (Q_former * eta_dot[i - 1 + index] * eta_dot[i - 1 + index] - 1 / sim.thermo_beta) / Qi;
        eta_dot[i + index] += eta_dot_dot[i + index] * dt4;
        eta_dot[i + index] *= exp_factor;
        Q_former = Qi;
    }

    // Update the derivatives of eta for the last component
    eta_dot_dot[nchains - 1 + index] = (Qi * eta_dot[nchains - 2 + index] * eta_dot[nchains - 2 + index] - 1 / sim.thermo_beta) / Qi;
    eta_dot[nchains - 1 + index] += eta_dot_dot[nchains - 1 + index] * dt4;

    return scale;
}

/* -------------------------------- */

NoseHooverNpThermostat::NoseHooverNpThermostat(Simulation& _sim, bool normal_modes, int _nchains) :
    NoseHooverThermostat(_sim, normal_modes, _nchains) {
    Q1 = NDIM * Qi;
    int eta_length = nchains * sim.natoms;
    eta = std::vector<double>(eta_length, 0.0);
    eta_dot = std::vector<double>(eta_length, 0.0);
    eta_dot_dot = std::vector<double>(eta_length, 0.0);

    required_energy = NDIM / sim.thermo_beta;
}

void NoseHooverNpThermostat::momentaUpdate() {
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        // Calculates the current energy in the system
        double current_energy = 0.0;
        for (int axis = 0; axis < NDIM; ++axis) {
            double momentum_for_calc = coupling->getMomentumForCalc(ptcl_idx, axis);
            current_energy += momentum_for_calc * momentum_for_calc;
        }
        current_energy /= sim.mass;

        // Obtain the scaling factor
        double scale = singleChainStep(current_energy, ptcl_idx * nchains);

        // Rescale the momenta
        for (int axis = 0; axis < NDIM; ++axis) {
            double momentum_for_calc = coupling->getMomentumForCalc(ptcl_idx, axis);
            double& momentum_for_update = coupling->getMomentumForUpdate(ptcl_idx, axis);
            momentum_for_update = momentum_for_calc * scale;
        }
    }
}

double NoseHooverNpThermostat::getAdditionToH() {
    double additionToH = 0;
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        additionToH += singleChainGetAdditionToH(NDIM, ptcl_idx * nchains);
    }
    return additionToH;
}

/* -------------------------------- */

NoseHooverNpDimThermostat::NoseHooverNpDimThermostat(Simulation& _sim, bool normal_modes, int _nchains) :
    NoseHooverThermostat(_sim, normal_modes, _nchains) {
    Q1 = Qi;

    int eta_length = nchains * sim.natoms * NDIM;
    eta = std::vector<double>(eta_length, 0.0);
    eta_dot = std::vector<double>(eta_length, 0.0);
    eta_dot_dot = std::vector<double>(eta_length, 0.0);

    required_energy = 1 / sim.thermo_beta;
}

void NoseHooverNpDimThermostat::momentaUpdate() {
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double momentum_for_calc = coupling->getMomentumForCalc(ptcl_idx, axis);
            // Calculates the current energy in the system
            double current_energy = momentum_for_calc * momentum_for_calc / sim.mass;

            // Obtain the scaling factor
            double scale = singleChainStep(current_energy, (ptcl_idx * NDIM + axis) * nchains);

            // Rescale the momentum
            double& momentum_for_update = coupling->getMomentumForUpdate(ptcl_idx, axis);
            momentum_for_update = momentum_for_calc * scale;
        }
    }
}

double NoseHooverNpDimThermostat::getAdditionToH() {
    double additionToH = 0;
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            additionToH += singleChainGetAdditionToH(1, (ptcl_idx * NDIM + axis) * nchains);
        }
    }
    return additionToH;
}

/* -------------------------------- */
