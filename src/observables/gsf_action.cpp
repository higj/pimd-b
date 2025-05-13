#include "observables/gsf_action.h"
#include "potentials.h"
#include "units.h"
#include <ranges>

/**
 * @brief Constructor for the class handling observables associated with the GSF action.
 */
GSFActionObservable::GSFActionObservable(Params& param_obj, int _freq, const std::string& _out_unit, int this_bead,
                                         Potential& ext_potential, Potential& int_potential, dVec& coord) :
    Observable(param_obj, _freq, _out_unit, this_bead), ext_potential(ext_potential), int_potential(int_potential), coord(coord) {
    getVariant(param_obj.sys["natoms"], natoms);
    getVariant(param_obj.sim["nbeads"], nbeads);
    getVariant(param_obj.sim["pbc"], pbc);
    getVariant(param_obj.sys["size"], size);
    
    double temperature;
    getVariant(param_obj.sys["temperature"], temperature);
    beta = 1.0 / (Constants::kB * temperature);

    std::string interaction_potential_name = std::get<std::string>(param_obj.interaction_pot["name"]);
    int_pot_cutoff = (interaction_potential_name == "free") ? 0.0 : std::get<double>(param_obj.interaction_pot["cutoff"]);

    double mass;
    getVariant(param_obj.sys["mass"], mass);
    double omega_p = nbeads / (beta * Constants::hbar);
    spring_constant = mass * omega_p * omega_p;

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

    double total_potential = ext_potential.V(coord);

    dVec gradients(natoms);
    gradients = ext_potential.gradV(coord);

    if (int_pot_cutoff != 0.0) {
        for (int ptcl_one = 0; ptcl_one < natoms; ++ptcl_one) {
            for (int ptcl_two = ptcl_one + 1; ptcl_two < natoms; ++ptcl_two) {
                dVec diff = getSeparation(ptcl_one, ptcl_two, MINIM, pbc, coord, size);  // Vectorial distance

                if (const double distance = diff.norm(); distance < int_pot_cutoff || int_pot_cutoff < 0.0) {
                    total_potential += int_potential.V(diff);
                    gradients = gradients + int_potential.gradV(diff);
                }
            }
        }
    }

    double total_force_squared = 0.0;

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            total_force_squared += gradients(ptcl_idx, axis) * gradients(ptcl_idx, axis);
        }
    }

    // Ensure the spring constant is in Tuckerman's convention
    
#if IPI_CONVENTION
    double sp_constant = spring_constant / nbeads;
#else
    double sp_constant = spring_constant;
#endif

    double potential_term = total_potential / (3 * nbeads);
    double force_squared_term = total_force_squared / (9 * sp_constant * nbeads * nbeads);

    if (this_bead % 2 != 0) {
        // Odd
        quantities["w_gsf"] = (-1.0) * potential_term + alpha * force_squared_term;

        // Evaluate the potential energy estimator only for odd imaginary-time slices
        quantities["pot_gsf"] = Units::convertToUser("energy", out_unit, total_potential / (0.5 * nbeads));
    } else {
        // Even
        quantities["w_gsf"] = potential_term + (1 - alpha) * force_squared_term;
    }

    quantities["w_gsf"] *= (-1.0) * beta;
}
