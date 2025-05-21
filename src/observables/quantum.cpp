#include "observables/quantum.h"
#include "bosonic_exchange.h"
#include "potentials.h"
#include "units.h"
#include <ranges>

/**
 * @brief Energy observable class constructor.
 */
QuantumObservable::QuantumObservable(Params& param_obj, int _freq, const std::string& _out_unit, int this_bead, 
                                     dVec& prev_coord, dVec& coord, BosonicExchangeBase& bosonic_exchange,
                                     Potential& ext_potential, Potential& int_potential, dVec& physical_forces, int& md_step) :
    EnergyObservable(param_obj, _freq, _out_unit, this_bead, prev_coord, coord, bosonic_exchange),
    ext_potential(ext_potential), int_potential(int_potential), physical_forces(physical_forces), md_step(md_step) {
        
    external_potential_name = std::get<std::string>(param_obj.external_pot["name"]);
    interaction_potential_name = std::get<std::string>(param_obj.interaction_pot["name"]);
    if (external_potential_name == "free" && interaction_potential_name == "free") {
        initialize({ "kinetic" });
    } else if (external_potential_name == "free" || interaction_potential_name == "free") {
        initialize({ "kinetic", "potential", "virial" });
    } else {
        initialize({ "kinetic", "potential", "ext_pot", "int_pot", "virial" });
    }
    double temperature;
    getVariant(param_obj.sys["temperature"], temperature);
    beta = 1.0 / (Constants::kB * temperature);
}

void QuantumObservable::calculate() {
    calculateKinetic();
    calculatePotential();
}

/**
 * @brief Calculates the quantum kinetic energy of the system using the primitive kinetic energy estimator.
 * Works both for distinguishable particles and bosons.
 */
void QuantumObservable::calculateKinetic() {
    // First, add the constant factor of d*N*P/(2*beta) to the kinetic energy (per bead)
    quantities["kinetic"] = 0.5 * NDIM * natoms / beta;

    // Then, subtract the spring energies. In the case of bosons, the exterior
    // spring energy requires separate treatment.
    if (this_bead == 0 && bosonic) {
        quantities["kinetic"] += bosonic_exchange.primEstimator();
    } else {
        double spring_energy = classicalSpringEnergy();
#if IPI_CONVENTION
        spring_energy /= nbeads;
#endif

        quantities["kinetic"] -= spring_energy;
    }

    quantities["kinetic"] = Units::convertToUser("energy", out_unit, quantities["kinetic"]);
}

/**
 * @brief Calculates the quantum potential energy of the system, based on the potential energy estimator.
 * The potential energy is the sum of the external potential energy and the interaction potential energy
 * across all time-slices, divided by the number of beads. In addition, the method calculates the virial
 * kinetic energy of the system.
 */
void QuantumObservable::calculatePotential() {
    double potential = 0.0;  // Total potential energy
    double virial = 0.0;     // Virial kinetic energy
    double int_pot = 0.0;    // Potential energy due to interactions
    double ext_pot = 0.0;    // Potential energy due to external field

    if (external_potential_name != "free") {
        ext_pot = ext_potential.getV(coord, md_step);
        potential += ext_pot;

        for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                virial -= coord(ptcl_idx, axis) * physical_forces(ptcl_idx, axis);
            }
        }
    }

    if (external_potential_name != "free" && interaction_potential_name != "free") {
        ext_pot /= nbeads;
        int_pot /= nbeads;

        quantities["ext_pot"] = Units::convertToUser("energy", out_unit, ext_pot);
        quantities["int_pot"] = Units::convertToUser("energy", out_unit, int_pot);
    }

    if (external_potential_name != "free" || interaction_potential_name != "free") {
        potential /= nbeads;
        virial *= 0.5 / nbeads;

        quantities["potential"] = Units::convertToUser("energy", out_unit, potential);
        quantities["virial"] = Units::convertToUser("energy", out_unit, virial);
    }
}