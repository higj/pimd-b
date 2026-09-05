#include "observables/suzuki_chin_virial_energy_observable.h"
#include "observables/shared_computations.h"
#include "core/force_manager.h"

SuzukiChinVirialEnergyObservable::SuzukiChinVirialEnergyObservable(
    const std::shared_ptr<const VecArray>& coord,
    const std::shared_ptr<const ForceManager>& force_mgr,
    const BeadContext& bead_ctx,
    int out_freq,
    const std::string& out_unit
) : Observable("virial_gsf", out_freq, out_unit),
m_coord(coord), m_force_mgr(force_mgr), m_bead_ctx(bead_ctx) {
    initializeLabel("virial_gsf");
}

void SuzukiChinVirialEnergyObservable::calculate() {
    if (m_bead_ctx.this_bead % 2 != 0) {  // odd beads only
        const auto [_, __, virial] =
            SharedComputations::suzukiChinComponents(*m_coord, *m_force_mgr, m_bead_ctx, m_cache);
        quantities["virial_gsf"] = Units::convertToUser(
            "energy", 
            m_out_unit, 
            virial / m_bead_ctx.nbeads
        );
    }
    // even beads: stay at 0.0
}