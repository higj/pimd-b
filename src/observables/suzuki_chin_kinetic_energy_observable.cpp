#include "observables/suzuki_chin_kinetic_energy_observable.h"
#include "observables/shared_computations.h"
#include "core/force_manager.h"

SuzukiChinKineticEnergyObservable::SuzukiChinKineticEnergyObservable(
    const std::shared_ptr<const VecArray>& coord,
    const std::shared_ptr<const ForceManager>& force_mgr,
    const BeadContext& bead_ctx,
    const BoxContext& box_ctx,
    const SpringContext& spring_ctx,
    double alpha,
    int out_freq,
    const std::string& out_unit
) : Observable("kin_gsf", out_freq, out_unit),
m_coord(coord), m_force_mgr(force_mgr),
m_bead_ctx(bead_ctx), m_box_ctx(box_ctx), m_spring_ctx(spring_ctx), m_alpha(alpha) {
    initializeLabel("kin_gsf");
}

void SuzukiChinKineticEnergyObservable::calculate() {
    const double force_sq = SharedComputations::suzukiChinComponents(
        *m_coord, *m_force_mgr, m_bead_ctx, m_box_ctx, m_cache
    ).force_squared;

#if TAU_CONVENTION
    const double sp_constant = m_spring_ctx.spring_constant / m_bead_ctx.nbeads;
#else
    const double sp_constant = m_spring_ctx.spring_constant;
#endif

    const double second_prefactor = 1.0 / (9.0 * sp_constant * m_bead_ctx.nbeads * m_bead_ctx.nbeads);
    const double alpha_eff = (m_bead_ctx.this_bead % 2 != 0) ? m_alpha : (1.0 - m_alpha);
    const double kinetic = second_prefactor * alpha_eff * force_sq;

    quantities["kin_gsf"] = Units::convertToUser("energy", m_out_unit, kinetic);
}