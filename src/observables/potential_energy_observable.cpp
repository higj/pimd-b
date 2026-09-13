#include "observables/potential_energy_observable.h"
#include "observables/shared_computations.h"
#include "core/force_manager.h"

PotentialEnergyObservable::PotentialEnergyObservable(
    const std::shared_ptr<const VecArray>& coord,
    const std::shared_ptr<const ForceManager>& force_mgr,
    const BoxContext& box_ctx,
    const BeadContext& bead_ctx,
    int out_freq,
    const std::string& out_unit
) : Observable("potential", out_freq, out_unit),
m_coord(coord),
m_force_mgr(force_mgr),
m_box_ctx(box_ctx),
m_bead_ctx(bead_ctx),
m_is_ext_free(force_mgr->ext_potential->isFree()),
m_is_int_free(force_mgr->int_potential->isFree()) {
    initializeLabel("potential");
}

void PotentialEnergyObservable::calculate() {
    double potential = 0.0;

    if (!m_is_ext_free) {
        potential += SharedComputations::extPotentialRaw(*m_coord, *m_force_mgr, m_cache);
    }

    if (!m_is_int_free) {
        potential += SharedComputations::intPotentialRaw(*m_coord, *m_force_mgr, m_box_ctx, m_bead_ctx, m_cache);
    }

    quantities["potential"] = Units::convertToUser(
        "energy",
        m_out_unit, 
        potential / m_bead_ctx.nbeads
    );
}