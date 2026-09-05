#include "observables/suzuki_chin_even_pot_energy_observable.h"
#include "observables/shared_computations.h"
#include "core/force_manager.h"

SuzukiChinEvenPotEnergyObservable::SuzukiChinEvenPotEnergyObservable(
    const std::shared_ptr<const VecArray>& coord,
    const std::shared_ptr<const ForceManager>& force_mgr,
    const BeadContext& bead_ctx,
    int out_freq,
    const std::string& out_unit
) : Observable("even_pot_gsf", out_freq, out_unit),
m_coord(coord), m_force_mgr(force_mgr), m_bead_ctx(bead_ctx) {
    initializeLabel("even_pot_gsf");
}

void SuzukiChinEvenPotEnergyObservable::calculate() {
    if (m_bead_ctx.this_bead % 2 == 0) {  // even beads only
        const auto [total_pot, _, __] =
            SharedComputations::suzukiChinComponents(*m_coord, *m_force_mgr, m_bead_ctx, m_cache);
        quantities["even_pot_gsf"] = Units::convertToUser(
            "energy", 
            m_out_unit, 
            total_pot / (0.5 * m_bead_ctx.nbeads)
        );
    }
    // odd beads: stay at 0.0; MPI sum accumulates even beads only
}