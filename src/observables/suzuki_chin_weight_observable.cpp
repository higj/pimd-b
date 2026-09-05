#include "observables/suzuki_chin_weight_observable.h"
#include "observables/shared_computations.h"
#include "core/force_manager.h"

SuzukiChinWeightObservable::SuzukiChinWeightObservable(
    const std::shared_ptr<const VecArray>& coord,
    const std::shared_ptr<const ForceManager>& force_mgr,
    const BeadContext& bead_ctx,
    const BoxContext& box_ctx,
    const ThermalContext& thermal_ctx,
    const SpringContext& spring_ctx,
    double alpha,
    int out_freq,
    const std::string& out_unit
) : Observable("w_gsf", out_freq, out_unit),
    m_coord(coord), m_force_mgr(force_mgr),
    m_bead_ctx(bead_ctx), m_box_ctx(box_ctx), m_thermal_ctx(thermal_ctx),
    m_spring_ctx(spring_ctx), m_alpha(alpha)
{
    initializeLabel("w_gsf");
}

void SuzukiChinWeightObservable::calculate()
{
    const auto [total_pot, force_sq, _] =
        SharedComputations::suzukiChinComponents(*m_coord, *m_force_mgr, m_bead_ctx, m_box_ctx, m_cache);

#if TAU_CONVENTION
    const double sp_constant = m_spring_ctx.spring_constant / m_bead_ctx.nbeads;
#else
    const double sp_constant = m_spring_ctx.spring_constant;
#endif

    const double potential_term = total_pot / (3.0 * m_bead_ctx.nbeads);
    const double force_sq_term = force_sq / (9.0 * sp_constant * m_bead_ctx.nbeads * m_bead_ctx.nbeads);

    double w;
    if (m_bead_ctx.this_bead % 2 != 0)
    {
        w = -potential_term + m_alpha * force_sq_term; // odd
    }
    else
    {
        w = potential_term + (1.0 - m_alpha) * force_sq_term; // even
    }

    quantities["w_gsf"] = -m_thermal_ctx.beta * w;
}
