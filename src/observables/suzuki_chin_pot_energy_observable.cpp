#include "observables/suzuki_chin_pot_energy_observable.h"
#include "observables/shared_computations.h"
#include "core/force_manager.h"

SuzukiChinPotEnergyObservable::SuzukiChinPotEnergyObservable(
    const std::shared_ptr<const VecArray>& coord,
    const std::shared_ptr<const ForceManager>& force_mgr,
    const BeadContext& bead_ctx,
    const BoxContext& box_ctx,
    int out_freq,
    const std::string& out_unit
) : Observable("sc_pot", out_freq, out_unit),
    m_coord(coord), m_force_mgr(force_mgr), m_bead_ctx(bead_ctx), m_box_ctx(box_ctx)
{
    initializeLabel("pot_gsf");
}

void SuzukiChinPotEnergyObservable::calculate()
{
    if (m_bead_ctx.this_bead % 2 != 0)
    {
        // odd beads only
        const double total_pot = SharedComputations::suzukiChinComponents(
            *m_coord, *m_force_mgr, m_bead_ctx, m_box_ctx, m_cache
        ).total_potential;

        quantities["pot_gsf"] = Units::convertToUser(
            "energy",
            m_out_unit,
            total_pot / (0.5 * m_bead_ctx.nbeads)
        );
    }
    // even beads: quantities["pot_gsf"] stays at 0.0; MPI sum accumulates odd beads only
}
