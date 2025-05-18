#include <utility>

#include "core/force_field_manager.h"
#include "potentials.h"
#include "bosonic_exchange/bosonic_exchange_base.h"

ForceFieldManager::ForceFieldManager(const std::shared_ptr<const SimulationConfig>& config) : m_config(config) {
    m_ext_potential = initializePotential(m_config->ext_pot_name, m_config->ext_pot_params);
    m_int_potential = initializePotential(m_config->int_pot_name, m_config->int_pot_params);

    // If the interaction potential is set to "free", then the cutoff distance is meaningless
    cutoff = (m_config->int_pot_name == "free") ? 0.0 : std::get<double>(m_config->int_pot_params.at("cutoff"));

    // For cubic cells with PBC, the cutoff distance must be no greater than L/2 for consistency with
    // the minimum image convention (see 1.6.3 in Allen & Tildesley).
    if (m_config->pbc) {
        cutoff = std::min(cutoff, 0.5 * m_config->box_size);
    }
}

/**
 * Updates the physical forces acting on the particles. This includes both the forces
 * due external potentials and the interaction forces between the particles.
 *
 * @param state Object representing the current state of the system, including forces acting on particles.
 */
void ForceFieldManager::updatePhysicalForces(SystemState& state) const
{
    // Calculate the external forces acting on the particles
    state.physical_forces = (-1.0) * m_ext_potential->gradV(state.coord);

    if (cutoff == 0.0)
        return;
    
    for (int ptcl_one = 0; ptcl_one < m_config->natoms; ++ptcl_one) {
        for (int ptcl_two = ptcl_one + 1; ptcl_two < m_config->natoms; ++ptcl_two) {
            // Get the vector distance between the two particles.
            // Here "diff" contains just one vector of dimension NDIM.
            /// TODO: ADD MINIM IMAGE!!
            dVec diff = state.coord.getSeparation(ptcl_one, ptcl_two);

            // If the distance between the particles exceeds the cutoff length
            // then we assume the interaction is negligible and do not bother
            // calculating the force.
            // We use the convention that when cutoff < 0 then the interaction is
            // calculated for all distances.
            if (const double distance = diff.norm(); distance < cutoff || cutoff < 0.0) {
                dVec force_on_one = (-1.0) * m_int_potential->gradV(diff);

                for (int axis = 0; axis < NDIM; ++axis) {
                    state.physical_forces(ptcl_one, axis) += force_on_one(0, axis);
                    state.physical_forces(ptcl_two, axis) -= force_on_one(0, axis);
                }
            }
        }
    }
}

/**
 * Updates the spring forces array.
 *
 * @param state Object representing the current state of the system, including forces acting on particles.
 * @param exchange_state Object representing the state of the exchange algorithm.
 */
void ForceFieldManager::updateSpringForces(SystemState& state, const ExchangeState& exchange_state) const {
    if (exchange_state.is_bosonic_bead) {
        // If the simulation is bosonic and the current bead is either 1 or P, we calculate
        // the exterior spring forces in the appropriate bosonic class.
        updateBosonicSpringForces(state, exchange_state);
        return;
    }

    // If particles are distinguishable, or if the current bead is an interior bead,
    // the force is calculated based on the standard expression for distinguishable particles.
    updateDistinguishableSpringForces(state);
}

/**
 * Updates the spring forces exerted on the beads.
 * In the distinguishable case, the force is given by Eqn. (12.6.4) in Tuckerman (1st ed).
 *
 * @param state Object representing the current state of the system, including forces acting on particles.
 */
void ForceFieldManager::updateDistinguishableSpringForces(SystemState& state) const {
    for (int ptcl_idx = 0; ptcl_idx < m_config->natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double diff_prev = state.prev_coord(ptcl_idx, axis) - state.coord(ptcl_idx, axis);
            double diff_next = state.next_coord(ptcl_idx, axis) - state.coord(ptcl_idx, axis);

#if MINIM
            if (m_config->pbc) {
                applyMinimumImage(diff_prev, m_config->box_size);
                applyMinimumImage(diff_next, m_config->box_size);
            }
#endif

            state.spring_forces(ptcl_idx, axis) = m_config->spring_constant * (diff_prev + diff_next);
        }
    }
}
/**
 * Updates the spring forces exerted on the beads.
 * In the bosonic case, by default, the forces are evaluated using the algorithm
 * described in https://doi.org/10.1063/5.0173749. It is also possible to perform the
 * bosonic simulation using the original (inefficient) algorithm, that takes into
 * account all the N! permutations, by setting FACTORIAL_BOSONIC_ALGORITHM to true.
 *
 * @param state Object representing the current state of the system, including forces acting on particles.
 * @param exchange_state Object representing the state of the exchange algorithm.
 */
void ForceFieldManager::updateBosonicSpringForces(SystemState& state, const ExchangeState& exchange_state)
{
    exchange_state.bosonic_exchange->prepare();
    exchange_state.bosonic_exchange->exteriorSpringForce(state.spring_forces);
}

/**
 * Initializes the potential based on the input parameters.
 *
 * @param potential_name Name of the potential.
 * @param potential_options Physical parameters of the potential.
 * @return Pointer to the initialized potential.
 */
std::unique_ptr<Potential> ForceFieldManager::initializePotential(const std::string& potential_name, const VariantMap& potential_options) const {
    if (potential_name == "free") {
        return std::make_unique<Potential>();
    }

    if (potential_name == "harmonic") {
        double omega = std::get<double>(potential_options.at("omega"));
        return std::make_unique<HarmonicPotential>(m_config->mass, omega);
    }

    if (potential_name == "double_well") {
        double strength = std::get<double>(potential_options.at("strength"));
        double loc = std::get<double>(potential_options.at("location"));
        return std::make_unique<DoubleWellPotential>(m_config->mass, strength, loc);
    }

    if (potential_name == "dipole") {
        double strength = std::get<double>(potential_options.at("strength"));
        return std::make_unique<DipolePotential>(strength);
    }

    if (potential_name == "cosine") {
        double amplitude = std::get<double>(potential_options.at("amplitude"));
        double phase = std::get<double>(potential_options.at("phase"));
        return std::make_unique<CosinePotential>(amplitude, m_config->box_size, phase);
    }

    if (potential_name == "aziz") {
        return std::make_unique<AzizPotential>();
    }

    return std::make_unique<Potential>();
}
