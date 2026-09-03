#include "observables/int_pot_observable.h"
#include "observables/shared_computations.h"
#include "core/force_manager.h"

IntPotObservable::IntPotObservable(
    const std::shared_ptr<const VecArray>& coord,
    const std::shared_ptr<const ForceManager>& force_mgr,
    const BoxContext& box_ctx,
    const BeadContext& bead_ctx,
    int out_freq,
    const std::string& out_unit
) : Observable("int_pot", out_freq, out_unit),
m_coord(coord),
m_force_mgr(force_mgr),
m_box_ctx(box_ctx),
m_bead_ctx(bead_ctx) {
    initializeLabel("int_pot");
}

void IntPotObservable::calculate() {
    const double raw = SharedComputations::intPotentialRaw(
        *m_coord, *m_force_mgr, m_box_ctx, m_bead_ctx, m_cache);
    quantities["int_pot"] = Units::convertToUser(
        "energy", m_out_unit, raw / m_bead_ctx.nbeads);
}