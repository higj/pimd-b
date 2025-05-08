#include "observables/gsf_action.h"
#include "simulation.h"
#include "units.h"
#include <ranges>
#include "mpi.h"

/**
 * @brief Constructor for the class handling observables associated with the GSF action.
 */
GSFActionObservable::GSFActionObservable(const Simulation& _sim, int _freq, const std::string& _out_unit) :
    Observable(_sim, _freq, _out_unit) {
    initialize({ "w_gsf", "pot_gsf" });
}

/**
 * @brief Calculates the natural logarithm of the weight associated with the GSF action,
 * which is used for re-weighting the observables [See J. Chem. Phys. 135, 064104 (2011)].
 * Also calculates the potential energy estimator using the operator method (only at odd imaginary-time slices).
 * Any additional estimators used with the GSF action should be calculated here.
 */
void GSFActionObservable::calculate() {
    double alpha = 0.0;

    double total_potential = sim.ext_potential->V(sim.coord);

    dVec gradients(sim.natoms);
    gradients = sim.ext_potential->gradV(sim.coord);

    if (sim.int_pot_cutoff != 0.0) {
        for (int ptcl_one = 0; ptcl_one < sim.natoms; ++ptcl_one) {
            for (int ptcl_two = ptcl_one + 1; ptcl_two < sim.natoms; ++ptcl_two) {
                dVec diff = sim.getSeparation(ptcl_one, ptcl_two, MINIM);  // Vectorial distance

                if (const double distance = diff.norm(); distance < sim.int_pot_cutoff || sim.int_pot_cutoff < 0.0) {
                    total_potential += sim.int_potential->V(diff);
                    gradients = gradients + sim.int_potential->gradV(diff);
                }
            }
        }
    }

    double total_force_squared = 0.0;

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            total_force_squared += gradients(ptcl_idx, axis) * gradients(ptcl_idx, axis);
        }
    }

    // Ensure the spring constant is in Tuckerman's convention
    
#if IPI_CONVENTION
    double sp_constant = sim.spring_constant / sim.nbeads;
#else
    double sp_constant = sim.spring_constant;
#endif

    double potential_term = total_potential / (3 * sim.nbeads);
    double force_squared_term = total_force_squared / (9 * sp_constant * sim.nbeads * sim.nbeads);

    if (sim.this_bead % 2 != 0) {
        // Odd
        quantities["w_gsf"] = (-1.0) * potential_term + alpha * force_squared_term;

        // Evaluate the potential energy estimator only for odd imaginary-time slices
        quantities["pot_gsf"] = Units::convertToUser("energy", out_unit, total_potential / (0.5 * sim.nbeads));
    } else {
        // Even
        quantities["w_gsf"] = potential_term + (1 - alpha) * force_squared_term;
    }

    quantities["w_gsf"] *= (-1.0) * sim.beta;
}
