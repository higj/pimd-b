#include "thermostat.h"

#include <numbers>

#include "simulation.h"
#include "common.h"

Thermostat::Thermostat(Simulation& _sim) : sim(_sim) {}
void Thermostat::step(){}

LangevinThermostat::LangevinThermostat(Simulation& _sim) : Thermostat(_sim) {
    a = exp(-0.5 * sim.gamma * sim.dt);

#if IPI_CONVENTION
    b = sqrt((1 - a * a) * sim.mass * sim.nbeads / sim.beta);
#else
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

NoseHooverThermostat::NoseHooverThermostat(Simulation& _sim, int _nchains) : Thermostat(_sim) {
    nchains = _nchains;
#if IPI_CONVENTION
    Qi = sim.nbeads / (sim.beta * sim.gamma * sim.gamma);
#else
    Qi = 1 / (sim.beta * sim.gamma * sim.gamma);
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
    required_kinetic_energy = 0.5 * NDIM * sim.natoms * sim.nbeads / sim.beta;
#else
    required_kinetic_energy = 0.5 * NDIM * sim.natoms / sim.beta;
#endif
}

void NoseHooverThermostat::step() {
    // To do: avoid repitition, the next 7 lines also exist on observable.cpp
    double kinetic_energy = 0.0;
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            kinetic_energy += sim.momenta(ptcl_idx, axis) * sim.momenta(ptcl_idx, axis);
        }
    }
    kinetic_energy *= 0.5 / sim.mass;

    double scale = singleChainStep(kinetic_energy, 0);

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            sim.momenta(ptcl_idx, axis) *= scale;
        }
    }
}

double NoseHooverThermostat::singleChainStep(const double& kinetic_energy, const int& index) {
    double exp_factor;
    eta_dot_dot[index] = (kinetic_energy - required_kinetic_energy) / Q1;

    eta_dot[nchains - 1 + index] += eta_dot_dot[nchains - 1 + index] * dt4;
    for (int i = nchains-2; i >= 0; i--) {
        exp_factor = exp(-dt8*eta_dot[i + 1 + index]);
        eta_dot[i + index] *= exp_factor;
        eta_dot[i + index] += eta_dot_dot[i + index] * dt4;
        eta_dot[i + index] *= exp_factor;
    }

    double scale = exp(-dt2 * eta_dot[index]);
    eta_dot_dot[index] = (kinetic_energy * scale * scale - required_kinetic_energy) / Q1;

#if CALC_ETA
    for (int i = 0; i < nchains; i++)
      eta[i + index] += dt2 * eta_dot[i + index];
#endif

    eta_dot[index] *= exp_factor;
    eta_dot[index] += eta_dot_dot[index] * dt4;
    eta_dot[index] *= exp_factor;

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
#if IPI_CONVENTION
    eta_dot_dot[nchains - 1 + index] = (Qi * eta_dot[nchains - 2 + index] * eta_dot[nchains - 2 + index] - sim.nbeads / sim.beta) / Qi;
#else
    eta_dot_dot[nchains - 1 + index] = (Qi * eta_dot[nchains - 2 + index] * eta_dot[nchains - 2 + index] - 1 / sim.beta) / Qi;
#endif
    eta_dot[nchains - 1 + index] += eta_dot_dot[nchains - 1 + index] * dt4;

    return scale;
}

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
    required_kinetic_energy = 0.5 * NDIM * sim.nbeads/ sim.beta;
#else
    required_kinetic_energy = 0.5 * NDIM / sim.beta;
#endif
}

void NoseHooverNpThermostat::step() {
    
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        double kinetic_energy = 0.0;
        for (int axis = 0; axis < NDIM; ++axis) {
            kinetic_energy += sim.momenta(ptcl_idx, axis) * sim.momenta(ptcl_idx, axis);
        }
        kinetic_energy *= 0.5 / sim.mass;
        
        double scale = singleChainStep(kinetic_energy, ptcl_idx * sim.natoms);

        for (int axis = 0; axis < NDIM; ++axis) {
            sim.momenta(ptcl_idx, axis) *= scale;
        }
    }   
}

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
    required_kinetic_energy = 0.5 * sim.nbeads / sim.beta;
#else    
    required_kinetic_energy = 0.5 / sim.beta;
#endif
}

void NoseHooverNpDimThermostat::step() {
    
    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double kinetic_energy = 0.5 * sim.momenta(ptcl_idx, axis) * sim.momenta(ptcl_idx, axis) / sim.mass;
            double scale = singleChainStep(kinetic_energy, (ptcl_idx * NDIM + axis) * sim.natoms);
            sim.momenta(ptcl_idx, axis) *= scale;
        }
    }   
}
