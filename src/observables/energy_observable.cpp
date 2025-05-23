#include "observables/energy_observable.h"
#include "core/exchange_state.h"
#include "core/force_field_manager.h"
#include "bosonic_exchange/bosonic_exchange_base.h"
#include "ring_polymer_utils.h"

/**
 * @brief Energy observable class constructor.
 */
EnergyObservable::EnergyObservable(const EnergyObservableContext& obs_context, int out_freq, const std::string& out_unit) :
    Observable(out_freq, out_unit), m_context(obs_context) {
    if (obs_context.ext_pot_name == "free" && obs_context.int_pot_name == "free") {
        initialize({ "kinetic" });
    } else if (obs_context.ext_pot_name == "free" || obs_context.int_pot_name == "free") {
        initialize({ "kinetic", "potential", "virial" });
    } else {
        initialize({ "kinetic", "potential", "ext_pot", "int_pot", "virial" });
    }
}

void EnergyObservable::calculate() {
    calculateKinetic();
    calculatePotential();
}

/**
 * @brief Calculates the quantum kinetic energy of the system using the primitive kinetic energy estimator.
 * Works both for distinguishable particles and bosons.
 */
void EnergyObservable::calculateKinetic() {
    // First, add the constant factor of d*N*P/(2*beta) to the kinetic energy (per bead)
    quantities["kinetic"] = 0.5 * NDIM * m_context.natoms / m_context.beta;

    // Then, subtract the spring energies. In the case of bosons, the exterior
    // spring energy requires separate treatment.
    if (m_context.this_bead == 0 && m_context.bosonic) {
        quantities["kinetic"] += m_context.exchange_state->bosonic_exchange->primEstimator();
    } else {
        /// TODO: Think about best way to pass minimum_image and box_size. Presumably, we need to pass Box object here?
        double spring_energy = RingPolymerUtils::classicalSpringEnergy(
            *m_context.coord,
            *m_context.prev_coord, 
            m_context.spring_constant,
            m_context.pbc && MINIM, /// TODO: MINIM should become a parameter (mic_spring and mic_potential)
            m_context.box_size
        );
#if IPI_CONVENTION
        spring_energy /= m_context.nbeads;
#endif

        quantities["kinetic"] -= spring_energy;
    }

    quantities["kinetic"] = Units::convertToUser("energy", m_out_unit, quantities["kinetic"]);
}

/**
 * @brief Calculates the quantum potential energy of the system, based on the potential energy estimator.
 * The potential energy is the sum of the external potential energy and the interaction potential energy
 * across all time-slices, divided by the number of beads. In addition, the method calculates the virial
 * kinetic energy of the system.
 */
void EnergyObservable::calculatePotential() {
    double potential = 0.0;  // Total potential energy
    double virial = 0.0;     // Virial kinetic energy
    double int_pot = 0.0;    // Potential energy due to interactions
    double ext_pot = 0.0;    // Potential energy due to external field

    const auto& coord = *m_context.coord;

    if (m_context.ext_pot_name != "free") {
        ext_pot = m_context.force_mgr->m_ext_potential->V(coord);
        potential += ext_pot;

        dVec physical_forces(m_context.natoms);
        physical_forces = (-1.0) * m_context.force_mgr->m_ext_potential->gradV(coord);

        for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                virial -= coord(ptcl_idx, axis) * physical_forces(ptcl_idx, axis);
            }
        }
    }

    if (m_context.force_mgr->cutoff != 0.0) {
        for (int ptcl_one = 0; ptcl_one < m_context.natoms; ++ptcl_one) {
            for (int ptcl_two = ptcl_one + 1; ptcl_two < m_context.natoms; ++ptcl_two) {
                dVec diff = coord.getSeparation(ptcl_one, ptcl_two);  // Vectorial distance

                /// TODO: MINIM should become a parameter (mic_spring and mic_potential)
                if (m_context.pbc && MINIM) {
                    applyMinimumImage(diff, m_context.box_size);
                }

                if (const double distance = diff.norm(); distance < m_context.force_mgr->cutoff || m_context.force_mgr->cutoff < 0.0) {
                    dVec force_on_one = (-1.0) * m_context.force_mgr->m_int_potential->gradV(diff);

                    double int_pot_val = m_context.force_mgr->m_int_potential->V(diff);
                    potential += int_pot_val;
                    int_pot += int_pot_val;

                    for (int axis = 0; axis < NDIM; ++axis) {
                        virial -= coord(ptcl_one, axis) * force_on_one(0, axis);
                    }
                }
            }
        }
    }
    
    if (m_context.ext_pot_name != "free" && m_context.int_pot_name != "free") {
        ext_pot /= m_context.nbeads;
        int_pot /= m_context.nbeads;

        quantities["ext_pot"] = Units::convertToUser("energy", m_out_unit, ext_pot);
        quantities["int_pot"] = Units::convertToUser("energy", m_out_unit, int_pot);
    }

    if (m_context.ext_pot_name != "free" || m_context.int_pot_name != "free") {
        potential /= m_context.nbeads;
        virial *= 0.5 / m_context.nbeads;

        quantities["potential"] = Units::convertToUser("energy", m_out_unit, potential);
        quantities["virial"] = Units::convertToUser("energy", m_out_unit, virial);
    }
}