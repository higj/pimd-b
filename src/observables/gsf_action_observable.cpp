#include "observables/gsf_action_observable.h"
#include "core/force_manager.h"

GSFActionObservable::GSFActionObservable(
    const std::shared_ptr<const VecArray>& coord,
    const std::shared_ptr<const ForceManager>& force_mgr,
    const BeadContext& bead_ctx,
    const ThermalContext& thermal_ctx,
    const SpringContext& spring_ctx,
    int out_freq,
    const std::string& out_unit
) : Observable("gsf", out_freq, out_unit),
    m_coord(coord),
    m_force_mgr(force_mgr),
    m_bead_ctx(bead_ctx),
    m_thermal_ctx(thermal_ctx),
    m_spring_ctx(spring_ctx)
{
    initialize({"w_gsf", "pot_gsf"});
}

void GSFActionObservable::calculate()
{
    double alpha = 0.0;
    const auto& coord = *m_coord;
    double total_potential = m_force_mgr->ext_potential->V(coord);

    VecArray gradients(m_bead_ctx.natoms);
    gradients = m_force_mgr->ext_potential->gradV(coord);

    if (m_force_mgr->cutoff != 0.0)
    {
        for (int ptcl_one = 0; ptcl_one < m_bead_ctx.natoms; ++ptcl_one)
        {
            for (int ptcl_two = ptcl_one + 1; ptcl_two < m_bead_ctx.natoms; ++ptcl_two)
            {
                /// TODO: ADD MINIM IMAGE!!
                SingleVec diff = coord.getSeparationArray(ptcl_one, ptcl_two);

                if (const double distance = norm(diff); distance < m_force_mgr->cutoff || m_force_mgr->cutoff < 0.0)
                {
                    total_potential += m_force_mgr->int_potential->V(diff);
                    
                    SingleVec int_grad = m_force_mgr->int_potential->gradV(diff);
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        gradients(ptcl_one, axis) += int_grad[axis];
                        gradients(ptcl_two, axis) -= int_grad[axis];
                    }
                }
            }
        }
    }

    double total_force_squared = 0.0;

    for (int ptcl_idx = 0; ptcl_idx < m_bead_ctx.natoms; ++ptcl_idx)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            total_force_squared += gradients(ptcl_idx, axis) * gradients(ptcl_idx, axis);
        }
    }

    // Ensure the spring constant is in Tuckerman's convention
#if IPI_CONVENTION
    double sp_constant = m_spring_ctx.spring_constant / m_bead_ctx.nbeads;
#else
    double sp_constant = m_context.spring_constant;
#endif

    const double potential_term = total_potential / (3 * m_bead_ctx.nbeads);
    const double force_squared_term = total_force_squared / (9 * sp_constant * m_bead_ctx.nbeads * m_bead_ctx.nbeads);

    if (m_bead_ctx.this_bead % 2 != 0)
    {
        // Odd
        quantities["w_gsf"] = (-1.0) * potential_term + alpha * force_squared_term;

        // Evaluate the potential energy estimator only for odd imaginary-time slices
        quantities["pot_gsf"] = Units::convertToUser("energy", m_out_unit, total_potential / (0.5 * m_bead_ctx.nbeads));
    }
    else
    {
        // Even
        quantities["w_gsf"] = potential_term + (1 - alpha) * force_squared_term;
    }

    quantities["w_gsf"] *= (-1.0) * m_thermal_ctx.beta;
}
