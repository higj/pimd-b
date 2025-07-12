#pragma once

#include "observables/observable.h"
#include "common.h"
#include "contexts/bead_context.h"
#include "contexts/spring_context.h"
#include "contexts/thermal_context.h"

class ForceManager;

class GSFActionObservable final : public Observable
{
public:
    /**
     * @brief Constructor for the class handling observables associated with the GSF action.
     */
    GSFActionObservable(
        const std::shared_ptr<const dVec>& coord,
        const std::shared_ptr<const ForceManager>& force_mgr,
        const BeadContext& bead_ctx,
        const ThermalContext& thermal_ctx,
        const SpringContext& spring_ctx,
        int out_freq,
        const std::string& out_unit
    );

    /**
     * @brief Calculates the natural logarithm of the weight associated with the GSF action,
     * which is used for re-weighting the observables [See J. Chem. Phys. 135, 064104 (2011)].
     * Also calculates the potential energy estimator using the operator method (only at odd imaginary-time slices).
     * Any additional estimators used with the GSF action should be calculated here.
     */
    void calculate() override;

private:
    std::shared_ptr<const dVec> m_coord;
    std::shared_ptr<const ForceManager> m_force_mgr;

    BeadContext m_bead_ctx;
    ThermalContext m_thermal_ctx;
    SpringContext m_spring_ctx;
};
