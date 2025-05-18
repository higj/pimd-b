#pragma once

#include "observables/observable.h"
#include "contexts/observables/bosonic_observable_context.h"

class BosonicObservable final : public Observable {
public:
    BosonicObservable(const BosonicObservableContext& obs_context, int out_freq, const std::string& out_unit);

    void calculate() override;
private:
    BosonicObservableContext m_context;
};