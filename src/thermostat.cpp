#include "thermostat.h"
#include <numbers>

#include "simulation.h"
#include "common.h"
#include "normal_modes.h"

Thermostat::Thermostat(Simulation& _sim) : sim(_sim) {}
void Thermostat::step(){}

/* -------------------------------- */

LangevinThermostat::LangevinThermostat(Simulation& _sim) : Thermostat(_sim) {
    // CR: this->a
    a = exp(-0.5 * sim.gamma * sim.dt);
#if IPI_CONVENTION
    // CR: this->b
    b = sqrt((1 - a * a) * sim.mass * sim.nbeads / sim.beta);
#else
    // CR: this->b
    b = sqrt((1 - a * a) * sim.mass / sim.beta);
#endif

}

void LangevinThermostat::step() {
    //std::normal_distribution<double> normal; // At the moment we use the Marsaglia generator

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double noise = sim.mars_gen->gaussian(); // LAMMPS pimd/langevin random numbers

            // Perturb the momenta with a Langevin thermostat
            sim.momenta(ptcl_idx, axis) = a * sim.momenta(ptcl_idx, axis) + b * noise;
        }
    }
}

/* -------------------------------- */

LangevinThermostatNM::LangevinThermostatNM(Simulation& _sim) : LangevinThermostat(_sim) {}

void LangevinThermostatNM::step() {
    // MPI communcation
    sim.normal_modes->shareData();
    MPI_Barrier(MPI_COMM_WORLD);

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double noise = sim.mars_gen->gaussian(); // LAMMPS pimd/langevin random numbers
            const int glob_idx = sim.normal_modes->globIndexAtom(axis, ptcl_idx);
            // Perturb the normal modes momenta with a Langevin thermostat
            // CR: try to separate concerns between thermostat, normal modes
//            arr_momenta_nm = cartesian_to_nm(momenta);
//            arr_momenta_nm *= whatever;
//            arr_momenta = nm_to_cartesian_to_nm(momenta)
            sim.normal_modes->arr_momenta_nm[glob_idx + sim.this_bead] = a * sim.normal_modes->momentumCarToNM(glob_idx) + b * noise;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    sim.normal_modes->updateCartesianMomenta();
}

/* -------------------------------- */

NoseHooverThermostat::NoseHooverThermostat(Simulation& _sim, int _nchains) : Thermostat(_sim) {
    nchains = _nchains;
#if IPI_CONVENTION
    Qi = Constants::hbar * Constants::hbar * sim.beta;
#else
    Qi = Constants::hbar * Constants::hbar * sim.beta / sim.nbeads;
#endif
    Q1 = NDIM * sim.natoms * Qi;
#if CALC_ETA
    eta = std::vector<double>(nchains);
#endif
    eta_dot = std::vector<double>(nchains);
    eta_dot_dot = std::vector<double>(nchains);
    for (int i = 0; i < nchains; i++) {
#if CALC_ETA
        eta[i] = 0.0;
#endif
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

void NoseHooverThermostat::step() {
    // Calculates the current energy in the system
    double current_energy = 0.0;
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            current_energy += sim.momenta(ptcl_idx, axis) * sim.momenta(ptcl_idx, axis);
        }
    }
    current_energy *= 1 / sim.mass;

    // Obtain the scaling factor
    double scale = singleChainStep(current_energy, 0);

    // Rescale the momenta
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            sim.momenta(ptcl_idx, axis) *= scale;
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

#if CALC_ETA
    for (int i = 0; i < nchains; i++)
      eta[i + index] += dt2 * eta_dot[i + index];
#endif

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

NoseHooverThermostatNM::NoseHooverThermostatNM(Simulation& _sim, int _nchains) : NoseHooverThermostat(_sim, _nchains) {}
void NoseHooverThermostatNM::step() {
    // MPI communcation
    sim.normal_modes->shareData();
    MPI_Barrier(MPI_COMM_WORLD);
    
    iVec glob_idxs(sim.natoms);
    dVec momenta_nm(sim.natoms);

    // Calculates the current energy in the system
    double current_energy = 0.0;
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            glob_idxs(ptcl_idx, axis) = sim.normal_modes->globIndexAtom(axis, ptcl_idx);
            momenta_nm(ptcl_idx, axis) = sim.normal_modes->momentumCarToNM(glob_idxs(ptcl_idx, axis));
            current_energy += momenta_nm(ptcl_idx, axis) * momenta_nm(ptcl_idx, axis);
        }
    }
    current_energy *= 1 / sim.mass;

    // Obtain the scaling factor
    double scale = singleChainStep(current_energy, 0);

    // Rescale the momenta in normal modes
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            sim.normal_modes->arr_momenta_nm[glob_idxs(ptcl_idx, axis) + sim.this_bead] = scale * momenta_nm(ptcl_idx, axis);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    sim.normal_modes->updateCartesianMomenta();
}

/* -------------------------------- */

NoseHooverNpThermostat::NoseHooverNpThermostat(Simulation& _sim, int _nchains) : NoseHooverThermostat(_sim, _nchains) {
    Q1 = NDIM * Qi;
#if CALC_ETA
    eta = std::vector<double>(nchains * sim.natoms);
#endif
    eta_dot = std::vector<double>(nchains * sim.natoms);
    eta_dot_dot = std::vector<double>(nchains * sim.natoms);
    for (int i = 0; i < nchains * sim.natoms; i++) {
#if CALC_ETA
        eta[i] = 0.0;
#endif
        eta_dot[i] = 0.0;
        eta_dot_dot[i] = 0.0;
    }
#if IPI_CONVENTION
    required_energy = NDIM * sim.nbeads/ sim.beta;
#else
    required_energy = NDIM / sim.beta;
#endif
}

void NoseHooverNpThermostat::step() {
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        // Calculates the current energy in the system
        double current_energy = 0.0;
        for (int axis = 0; axis < NDIM; ++axis) {
            current_energy += sim.momenta(ptcl_idx, axis) * sim.momenta(ptcl_idx, axis);
        }
        current_energy *= 1 / sim.mass;
        
        // Obtain the scaling factor
        double scale = singleChainStep(current_energy, ptcl_idx * nchains);

        // Rescale the momenta
        for (int axis = 0; axis < NDIM; ++axis) {
            sim.momenta(ptcl_idx, axis) *= scale;
        }
    }   
}

/* -------------------------------- */

NoseHooverNpThermostatNM::NoseHooverNpThermostatNM(Simulation& _sim, int _nchains) : NoseHooverNpThermostat(_sim, _nchains) {}

void NoseHooverNpThermostatNM::step() {
    // MPI communcation
    sim.normal_modes->shareData();
    MPI_Barrier(MPI_COMM_WORLD);

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        std::vector<int> glob_idxs(NDIM);
        std::vector<double> momenta_nm(NDIM);

        // Calculates the current energy in the system
        double current_energy = 0.0;
        for (int axis = 0; axis < NDIM; ++axis) {
            glob_idxs[axis] = sim.normal_modes->globIndexAtom(axis, ptcl_idx);
            momenta_nm[axis] = sim.normal_modes->momentumCarToNM(glob_idxs[axis]);
            current_energy += momenta_nm[axis] * momenta_nm[axis];
        }
        current_energy *= 1 / sim.mass;

        // Obtain the scaling factor
        double scale = singleChainStep(current_energy, ptcl_idx * nchains);

        // Rescale the momenta in normal modes
        for (int axis = 0; axis < NDIM; ++axis) {
            sim.normal_modes->arr_momenta_nm[glob_idxs[axis] + sim.this_bead] = scale * momenta_nm[axis];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    sim.normal_modes->updateCartesianMomenta();
}
/* -------------------------------- */

NoseHooverNpDimThermostat::NoseHooverNpDimThermostat(Simulation& _sim, int _nchains) : NoseHooverThermostat(_sim, _nchains) {
    Q1 = Qi;
#if CALC_ETA
    eta = std::vector<double>(nchains * sim.natoms * NDIM);
#endif
    eta_dot = std::vector<double>(nchains * sim.natoms * NDIM);
    eta_dot_dot = std::vector<double>(nchains * sim.natoms * NDIM);
    for (int i = 0; i < nchains * sim.natoms * NDIM; i++) {
#if CALC_ETA
        eta[i] = 0.0;
#endif
        eta_dot[i] = 0.0;
        eta_dot_dot[i] = 0.0;
    }
#if IPI_CONVENTION
    required_energy = sim.nbeads / sim.beta;
#else    
    required_energy = 1 / sim.beta;
#endif
}

void NoseHooverNpDimThermostat::step() {
    
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            // Calculates the current energy in the system
            double current_energy = sim.momenta(ptcl_idx, axis) * sim.momenta(ptcl_idx, axis) / sim.mass;

            // Obtain the scaling factor
            double scale = singleChainStep(current_energy, (ptcl_idx * NDIM + axis) * nchains);

            // Rescale the momentum
            sim.momenta(ptcl_idx, axis) *= scale;
        }
    }   
}

/* -------------------------------- */

NoseHooverNpDimThermostatNM::NoseHooverNpDimThermostatNM(Simulation& _sim, int _nchains) : NoseHooverNpDimThermostat(_sim, _nchains) {}

void NoseHooverNpDimThermostatNM::step() {

    // MPI communcation
    sim.normal_modes->shareData();
    MPI_Barrier(MPI_COMM_WORLD);
    
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            const int glob_idx = sim.normal_modes->globIndexAtom(axis, ptcl_idx);
            double momentum_nm = sim.normal_modes->momentumCarToNM(glob_idx);

            // Calculates the current energy in the system
            double current_energy = momentum_nm * momentum_nm / sim.mass;
 
            // Obtain the scaling factor
            double scale = singleChainStep(current_energy, (ptcl_idx * NDIM + axis) * nchains);
 
            // Rescale the momentum in normal modes
            sim.normal_modes->arr_momenta_nm[glob_idx + sim.this_bead] = scale * momentum_nm;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    sim.normal_modes->updateCartesianMomenta();
}
