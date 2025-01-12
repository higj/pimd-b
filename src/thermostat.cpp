#include "thermostat.h"
#include <numbers>
#include "simulation.h"
#include "common.h"
#include "normal_modes.h"
#include "thermostat_coupling.h"

Thermostat::Thermostat(Simulation& _sim, bool normal_modes) : sim(_sim) {
    // Choose coupling (Cartesian coords or normal modes of distinguishable ring polymers)
    if (normal_modes) {
        coupling = std::make_unique<NMCoupling>(_sim);
    }
    else {
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
void Thermostat::momentaUpdate(){}

double Thermostat::getAdditionToH() {
    return 0;
}

/* -------------------------------- */

LangevinThermostat::LangevinThermostat(Simulation& _sim, bool normal_modes) : Thermostat(_sim, normal_modes) {
    // CR: this->a
    friction_coefficient = exp(-0.5 * sim.gamma * sim.dt);
#if IPI_CONVENTION
    // CR: this->b
    noise_coefficient = sqrt((1 - friction_coefficient * friction_coefficient) * sim.mass * sim.nbeads / sim.beta);
#else
    // CR: this->b
    noise_coefficient = sqrt((1 - friction_coefficient * friction_coefficient) * sim.mass / sim.beta);
#endif

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

NoseHooverThermostat::NoseHooverThermostat(Simulation& _sim, bool normal_modes, int _nchains) : Thermostat(_sim, normal_modes) {
    nchains = _nchains;
#if IPI_CONVENTION
    Qi = Constants::hbar * Constants::hbar * sim.beta;
#else
    Qi = Constants::hbar * Constants::hbar * sim.beta / sim.nbeads;
#endif
    Q1 = NDIM * sim.natoms * Qi;
    eta = std::vector<double>(nchains);
    eta_dot = std::vector<double>(nchains);
    eta_dot_dot = std::vector<double>(nchains);
    for (int i = 0; i < nchains; i++) {
        eta[i] = 0.0;
        eta_dot[i] = 0.0;
        eta_dot_dot[i] = 0.0;
    }
    dt2 = 0.5 * sim.dt;
    dt4 = 0.25 * sim.dt;
    dt8 = 0.125 * sim.dt;
#if IPI_CONVENTION
    required_energy = NDIM * sim.natoms * sim.nbeads / sim.beta;
#else
    required_energy = NDIM * sim.natoms / sim.beta;
#endif
}

double NoseHooverThermostat::getAdditionToH() {
    return singleChainGetAdditionToH(NDIM * sim.natoms, 0);    
}

double NoseHooverThermostat::singleChainGetAdditionToH(const int& expected_energy, const int& index) {
    double additionToH = 0.5 * Q1 * eta_dot[index] * eta_dot[index];
#if IPI_CONVENTION
    additionToH += expected_energy * eta[index] * sim.nbeads / sim.beta;
#else
    additionToH += expected_energy * eta[index] / sim.beta;
#endif

    for (int i = 1; i < nchains; i++) {
        additionToH += 0.5 * Qi * eta_dot[index + i] * eta_dot[index + i];
#if IPI_CONVENTION
        additionToH += eta[index + i] * sim.nbeads / sim.beta;
#else
        additionToH += eta[index + i] / sim.beta;
#endif
    }

    return additionToH;
}

void NoseHooverThermostat::momentaUpdate() {
    // Calculates the current energy in the system
    double current_energy = 0.0;
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double momentum_for_calc = coupling->getMomentumForCalc(ptcl_idx, axis);
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
 * Mark E Tuckerman et al 2006 J. Phys. A: Math. Gen. 39 5629
 * and its implementation in LAMMPS
 * 
 * @param current_energy The current energy
 * @param index          The index of the first component in the chain
 * @return               The scaling factor 
 */
double NoseHooverThermostat::singleChainStep(const double& current_energy, const int& index) {
    double exp_factor;

    // Update the second derivative of eta for the first component
    eta_dot_dot[index] = (current_energy - required_energy) / Q1;

    // Update the first derivative of eta for the last component
    eta_dot[nchains - 1 + index] += eta_dot_dot[nchains - 1 + index] * dt4;

    // Update the first derivative of eta for all others components
    for (int i = nchains-2; i >= 0; i--) {
        exp_factor = exp(-dt8*eta_dot[i + 1 + index]);
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
#if IPI_CONVENTION
        eta_dot_dot[i + index] = (Q_former * eta_dot[i - 1 + index] * eta_dot[i - 1 + index] - sim.nbeads / sim.beta) / Qi;
#else
        eta_dot_dot[i + index] = (Q_former * eta_dot[i - 1 + index] * eta_dot[i - 1 + index] - 1 / sim.beta) / Qi;
#endif
	eta_dot[i + index] += eta_dot_dot[i + index] * dt4;
        eta_dot[i + index] *= exp_factor;
        Q_former = Qi;
    }

    // Update the derivatives of eta for the last component
#if IPI_CONVENTION
    eta_dot_dot[nchains - 1 + index] = (Qi * eta_dot[nchains - 2 + index] * eta_dot[nchains - 2 + index] - sim.nbeads / sim.beta) / Qi;
#else
    eta_dot_dot[nchains - 1 + index] = (Qi * eta_dot[nchains - 2 + index] * eta_dot[nchains - 2 + index] - 1 / sim.beta) / Qi;
#endif
    eta_dot[nchains - 1 + index] += eta_dot_dot[nchains - 1 + index] * dt4;

    return scale;
}

/* -------------------------------- */

NoseHooverNpThermostat::NoseHooverNpThermostat(Simulation& _sim, bool normal_modes, int _nchains) : NoseHooverThermostat(_sim, normal_modes, _nchains) {
    Q1 = NDIM * Qi;
    eta = std::vector<double>(nchains * sim.natoms);
    eta_dot = std::vector<double>(nchains * sim.natoms);
    eta_dot_dot = std::vector<double>(nchains * sim.natoms);
    for (int i = 0; i < nchains * sim.natoms; i++) {
        eta[i] = 0.0;
        eta_dot[i] = 0.0;
        eta_dot_dot[i] = 0.0;
    }
#if IPI_CONVENTION
    required_energy = NDIM * sim.nbeads/ sim.beta;
#else
    required_energy = NDIM / sim.beta;
#endif
}

void NoseHooverNpThermostat::momentaUpdate() {
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        // Calculates the current energy in the system
        double current_energy = 0.0;
        for (int axis = 0; axis < NDIM; ++axis) {
            double momentum_for_calc = coupling->getMomentumForCalc(ptcl_idx, axis);
            current_energy += momentum_for_calc * momentum_for_calc;
        }
        current_energy *= 1 / sim.mass;
        
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

NoseHooverNpDimThermostat::NoseHooverNpDimThermostat(Simulation& _sim, bool normal_modes, int _nchains) : NoseHooverThermostat(_sim, normal_modes, _nchains) {
    Q1 = Qi;
    eta = std::vector<double>(nchains * sim.natoms * NDIM);
    eta_dot = std::vector<double>(nchains * sim.natoms * NDIM);
    eta_dot_dot = std::vector<double>(nchains * sim.natoms * NDIM);
    for (int i = 0; i < nchains * sim.natoms * NDIM; i++) {
        eta[i] = 0.0;
        eta_dot[i] = 0.0;
        eta_dot_dot[i] = 0.0;
    }
#if IPI_CONVENTION
    required_energy = sim.nbeads / sim.beta;
#else    
    required_energy = 1 / sim.beta;
#endif
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
