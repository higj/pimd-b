#include "observables/energy.h"
#include "params.h"
#include <numbers>

/**
 * @brief Classical observable class constructor.
 */
EnergyObservable::EnergyObservable(Params& param_obj, int _freq, const std::string& _out_unit, 
                                   int this_bead, dVec& prev_coord, dVec& coord, BosonicExchangeBase& bosonic_exchange) :
    Observable(param_obj, _freq, _out_unit, this_bead), prev_coord(prev_coord), coord(coord), bosonic_exchange(bosonic_exchange) {
    getVariant(param_obj.sys["natoms"], natoms);
    getVariant(param_obj.sim["nbeads"], nbeads);
    getVariant(param_obj.sim["bosonic"], bosonic);
    getVariant(param_obj.sim["pbc"], pbc);
    getVariant(param_obj.sys["mass"], mass);
    getVariant(param_obj.sys["size"], size);
    double temperature;
    getVariant(param_obj.sys["temperature"], temperature);
#if IPI_CONVENTION
    // i-Pi convention [J. Chem. Phys. 133, 124104 (2010)]
    double omega_p = nbeads * Constants::kB * temperature / Constants::hbar;
#else
    // Tuckerman convention
    double omega_p = sqrt(nbeads) * Constants::kB * temperature / Constants::hbar;
#endif
    spring_constant = mass * omega_p * omega_p;
}

void EnergyObservable::calculate() {}

/**
 * Calculates the spring energy contribution of the current and the previous time-slice,
 * provided they are classical, i.e., are not affected by bosonic exchange.
 * If the simulation is bosonic, the function is callable only for the interior connections.
 *
 * @return Classical spring energy contribution of the current and the previous time-slice.
 */
double EnergyObservable::classicalSpringEnergy() const {
    assert(!bosonic || (bosonic && this_bead != 0));

    double interior_spring_energy = 0.0;

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double diff = prev_coord(ptcl_idx, axis) - coord(ptcl_idx, axis);

#if MINIM
            if (pbc) {
                applyMinimumImage(diff, size);
            }
#endif

            interior_spring_energy += diff * diff;
        }
    }

    interior_spring_energy *= 0.5 * spring_constant;
    return interior_spring_energy;
}