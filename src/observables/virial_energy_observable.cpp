#include "observables/virial_energy_observable.h"
#include "observables/shared_computations.h"

VirialEnergyObservable::VirialEnergyObservable(
    const std::shared_ptr<const VecArray>& coord,
    const std::shared_ptr<const VecArray>& physical_forces,
    const BeadContext& bead_ctx,
    int out_freq,
    const std::string& out_unit
) : Observable("virial", out_freq, out_unit),
m_coord(coord),
m_forces(physical_forces),
m_bead_ctx(bead_ctx) {
    initializeLabel("virial");
}

void VirialEnergyObservable::calculate() {
    const double raw = SharedComputations::virialRaw(*m_coord, *m_forces, m_bead_ctx, m_cache);
    quantities["virial"] = Units::convertToUser(
        "energy", 
        m_out_unit, 
        0.5 * raw / m_bead_ctx.nbeads
    );
}