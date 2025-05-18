#pragma once

#include "observables/observable.h"
#include "contexts/observables/energy_observable_context.h"

class EnergyObservable final : public Observable {
public:
    EnergyObservable(const EnergyObservableContext& obs_context, int out_freq, const std::string& out_unit);

    void calculate() override;

private:
    EnergyObservableContext m_context;

    void calculateKinetic();
    void calculatePotential();
};