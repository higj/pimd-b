#include "thermostats/nose_hoover_thermostat.h"
#include "thermostats/thermostat_coupling.h"
#include "common.h"

// A class for a single Nose-Hoover chain coupled to all degrees of freedom
NoseHooverThermostat::NoseHooverThermostat(const ThermostatContext& context, const NoseHooverThermostatContext& nh_context) : Thermostat(
    context), m_nh_context(nh_context) {
    // Q1 and Qi (where i>1) are the effective masses of the extended variables, which determine
    // the time scale on which the heat-bath variables evolve.
    // The choice is based on Sec. 2.5 in Martyna, G. J. et al. (1996), Mol. Phys., 87(5), pp. 1117-1157.
    // For the fluctuation frequency, we choose the spring frequency of the ring polymers.
    Qi = Constants::hbar * Constants::hbar * context.beta / context.nbeads;

    // Q1/Qi should be equal to the number of degrees of freedom (NDOF).
    // In the absence of constraints, NDOF=NDIM*natoms.
    Q1 = NDIM * context.natoms * Qi;

    eta = std::vector(m_nh_context.nchains, 0.0);
    eta_dot = std::vector(m_nh_context.nchains, 0.0);
    eta_dot_dot = std::vector(m_nh_context.nchains, 0.0);

    dt2 = 0.5 * context.dt;
    dt4 = 0.25 * context.dt;
    dt8 = 0.125 * context.dt;

    required_energy = NDIM * context.natoms / context.thermo_beta;
}

/**
 * Returns the contribution of NHC to the conserved quantity.
 * Corresponds to the difference between the conserved quantity
 * and the classical Hamiltonian of the ring polymers.
*/
 
double NoseHooverThermostat::getAdditionToH() {
    return singleChainGetAdditionToH(NDIM * m_context.natoms, 0);
}

/**
 * Calculate the contribution of NHC to the conserved quantity.
 *
 * @param ndof  Number of degrees of freedom
 * @param index Index of the first component in the chain
 * @return      Energy contribution to the conserved quantity
 */
double NoseHooverThermostat::singleChainGetAdditionToH(const int& ndof, const int& index) {
    // First, we add the kinetic energy associated with the thermostat.
    // Note that the momentum of the thermostat is given by Qi*eta_dot.
    double addition_to_h = 0.5 * Q1 * eta_dot[index] * eta_dot[index];

    // Then, we add the other term that is proportional to kT*eta
    addition_to_h += ndof * eta[index] / m_context.thermo_beta;

    // Similar calculation, but for the rest of the thermostats
    for (int i = 1; i < m_nh_context.nchains; i++) {
        addition_to_h += 0.5 * Qi * eta_dot[index + i] * eta_dot[index + i];
        addition_to_h += eta[index + i] / m_context.thermo_beta;
    }

    return addition_to_h;
}

void NoseHooverThermostat::momentaUpdate() {
    // Calculates the current energy in the system
    double current_energy = 0.0;
    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            const double momentum_for_calc = coupling->getMomentumForCalc(ptcl_idx, axis);
            current_energy += momentum_for_calc * momentum_for_calc;
        }
    }
    current_energy /= m_context.mass;

    // Obtain the scaling factor
    const double scale = singleChainStep(current_energy, 0);

    // Rescale the momenta
    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
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
    eta_dot[m_nh_context.nchains - 1 + index] += eta_dot_dot[m_nh_context.nchains - 1 + index] * dt4;

    // Update the first derivative of eta for all others components
    for (int i = m_nh_context.nchains - 2; i >= 0; i--) {
        exp_factor = exp(-dt8 * eta_dot[i + 1 + index]);
        eta_dot[i + index] *= exp_factor;
        eta_dot[i + index] += eta_dot_dot[i + index] * dt4;
        eta_dot[i + index] *= exp_factor;
    }

    // Calculate the rescaling factor
    double scale = exp(-dt2 * eta_dot[index]);
    for (int i = 0; i < m_nh_context.nchains; i++)
        eta[i + index] += dt2 * eta_dot[i + index];

    // Update the derivatives of eta for the first component
    eta_dot_dot[index] = (current_energy * scale * scale - required_energy) / Q1;
    eta_dot[index] *= exp_factor;
    eta_dot[index] += eta_dot_dot[index] * dt4;
    eta_dot[index] *= exp_factor;

    // Update the derivatives of eta for all other components except the last
    double Q_former = Q1;
    for (int i = 1; i < m_nh_context.nchains - 1; i++) {
        exp_factor = exp(-dt8 * eta_dot[i + 1 + index]);
        eta_dot[i + index] *= exp_factor;
        eta_dot_dot[i + index] = (Q_former * eta_dot[i - 1 + index] * eta_dot[i - 1 + index] - 1 / m_context.thermo_beta) / Qi;
        eta_dot[i + index] += eta_dot_dot[i + index] * dt4;
        eta_dot[i + index] *= exp_factor;
        Q_former = Qi;
    }

    // Update the derivatives of eta for the last component
    eta_dot_dot[m_nh_context.nchains - 1 + index] = (Qi * eta_dot[m_nh_context.nchains - 2 + index] * eta_dot[m_nh_context.nchains - 2 + index] - 1 / m_context.thermo_beta) / Qi;
    eta_dot[m_nh_context.nchains - 1 + index] += eta_dot_dot[m_nh_context.nchains - 1 + index] * dt4;

    return scale;
}

/* -------------------------------- */

// A class for a unique Nose-Hoover chain coupled to each particle
NoseHooverNpThermostat::NoseHooverNpThermostat(const ThermostatContext& context, const NoseHooverThermostatContext& nh_context) :
    NoseHooverThermostat(context, nh_context) {
    Q1 = NDIM * Qi;
    const int eta_length = nh_context.nchains * context.natoms;
    eta = std::vector(eta_length, 0.0);
    eta_dot = std::vector(eta_length, 0.0);
    eta_dot_dot = std::vector(eta_length, 0.0);

    required_energy = NDIM / context.thermo_beta;
}

void NoseHooverNpThermostat::momentaUpdate() {
    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        // Calculates the current energy in the system
        double current_energy = 0.0;
        for (int axis = 0; axis < NDIM; ++axis) {
            double momentum_for_calc = coupling->getMomentumForCalc(ptcl_idx, axis);
            current_energy += momentum_for_calc * momentum_for_calc;
        }
        current_energy /= m_context.mass;

        // Obtain the scaling factor
        double scale = singleChainStep(current_energy, ptcl_idx * m_nh_context.nchains);

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
    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        additionToH += singleChainGetAdditionToH(NDIM, ptcl_idx * m_nh_context.nchains);
    }
    return additionToH;
}

/* -------------------------------- */

// A class for a unique Nose-Hoover chain coupled to each degree of freedom
NoseHooverNpDimThermostat::NoseHooverNpDimThermostat(const ThermostatContext& context, const NoseHooverThermostatContext& nh_context) :
    NoseHooverThermostat(context, nh_context) {
    Q1 = Qi;

    const int eta_length = nh_context.nchains * context.natoms * NDIM;
    eta = std::vector(eta_length, 0.0);
    eta_dot = std::vector(eta_length, 0.0);
    eta_dot_dot = std::vector(eta_length, 0.0);

    required_energy = 1 / context.thermo_beta;
}

void NoseHooverNpDimThermostat::momentaUpdate() {
    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double momentum_for_calc = coupling->getMomentumForCalc(ptcl_idx, axis);
            // Calculates the current energy in the system
            double current_energy = momentum_for_calc * momentum_for_calc / m_context.mass;

            // Obtain the scaling factor
            double scale = singleChainStep(current_energy, (ptcl_idx * NDIM + axis) * m_nh_context.nchains);

            // Rescale the momentum
            double& momentum_for_update = coupling->getMomentumForUpdate(ptcl_idx, axis);
            momentum_for_update = momentum_for_calc * scale;
        }
    }
}

double NoseHooverNpDimThermostat::getAdditionToH() {
    double addition_to_h = 0;
    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            addition_to_h += singleChainGetAdditionToH(1, (ptcl_idx * NDIM + axis) * m_nh_context.nchains);
        }
    }
    return addition_to_h;
}

/* -------------------------------- */
