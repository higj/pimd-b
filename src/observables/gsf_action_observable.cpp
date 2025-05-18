#include "observables/gsf_action_observable.h"
#include "core/force_field_manager.h"

GSFActionObservable::GSFActionObservable(const GSFActionObservableContext& obs_context, int out_freq, const std::string& out_unit) :
    Observable(out_freq, out_unit), m_context(obs_context) {
    initialize({ "w_gsf", "pot_gsf" });
}

void GSFActionObservable::calculate() {
    double alpha = 0.0;
    const auto& coord = *m_context.coord;
    double total_potential = m_context.force_mgr->m_ext_potential->V(coord);

    dVec gradients(m_context.natoms);
    gradients = m_context.force_mgr->m_ext_potential->gradV(coord);

    if (m_context.force_mgr->cutoff != 0.0) {
        for (int ptcl_one = 0; ptcl_one < m_context.natoms; ++ptcl_one) {
            for (int ptcl_two = ptcl_one + 1; ptcl_two < m_context.natoms; ++ptcl_two) {
                /// TODO: ADD MINIM IMAGE!!
                dVec diff = coord.getSeparation(ptcl_one, ptcl_two);  // Vectorial distance

                if (const double distance = diff.norm(); distance < m_context.force_mgr->cutoff || m_context.force_mgr->cutoff < 0.0) {
                    total_potential += m_context.force_mgr->m_int_potential->V(diff);
                    gradients = gradients + m_context.force_mgr->m_int_potential->gradV(diff);
                }
            }
        }
    }

    double total_force_squared = 0.0;

    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            total_force_squared += gradients(ptcl_idx, axis) * gradients(ptcl_idx, axis);
        }
    }

    // Ensure the spring constant is in Tuckerman's convention
#if IPI_CONVENTION
    double sp_constant = m_context.spring_constant / m_context.nbeads;
#else
    double sp_constant = m_context.spring_constant;
#endif

    const double potential_term = total_potential / (3 * m_context.nbeads);
    const double force_squared_term = total_force_squared / (9 * sp_constant * m_context.nbeads * m_context.nbeads);

    if (m_context.this_bead % 2 != 0) {
        // Odd
        quantities["w_gsf"] = (-1.0) * potential_term + alpha * force_squared_term;

        // Evaluate the potential energy estimator only for odd imaginary-time slices
        quantities["pot_gsf"] = Units::convertToUser("energy", m_out_unit, total_potential / (0.5 * m_context.nbeads));
    } else {
        // Even
        quantities["w_gsf"] = potential_term + (1 - alpha) * force_squared_term;
    }

    quantities["w_gsf"] *= (-1.0) * m_context.beta;
}
