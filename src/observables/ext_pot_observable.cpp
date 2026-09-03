#include "observables/ext_pot_observable.h"
#include "observables/shared_computations.h"
#include "core/force_manager.h"

ExtPotObservable::ExtPotObservable(
    const std::shared_ptr<const VecArray>& coord,
    const std::shared_ptr<const ForceManager>& force_mgr,
    const BeadContext& bead_ctx,
    int out_freq,
    const std::string& out_unit
) : Observable("ext_pot", out_freq, out_unit),
m_coord(coord),
m_force_mgr(force_mgr),
m_bead_ctx(bead_ctx) {
    initializeLabel("ext_pot");
}

void ExtPotObservable::calculate() {
    const double raw = SharedComputations::extPotentialRaw(*m_coord, *m_force_mgr, m_cache);

    quantities["ext_pot"] = Units::convertToUser(
        "energy",
        m_out_unit, 
        raw / m_bead_ctx.nbeads
    );
}