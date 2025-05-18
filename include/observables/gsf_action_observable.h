#pragma once

#include "observables/observable.h"
#include "contexts/observables/gsf_action_observable_context.h"

class GSFActionObservable final : public Observable {
public:
    /**
     * @brief Constructor for the class handling observables associated with the GSF action.
     */
    GSFActionObservable(const GSFActionObservableContext& obs_context, int out_freq, const std::string& out_unit);

    /**
     * @brief Calculates the natural logarithm of the weight associated with the GSF action,
     * which is used for re-weighting the observables [See J. Chem. Phys. 135, 064104 (2011)].
     * Also calculates the potential energy estimator using the operator method (only at odd imaginary-time slices).
     * Any additional estimators used with the GSF action should be calculated here.
     */
    void calculate() override;
private:
    GSFActionObservableContext m_context;
};