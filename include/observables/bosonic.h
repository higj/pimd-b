#pragma once

#include "observables/observable.h"

class BosonicExchangeBase; // Forward declaration

/* -------------------------------- */

class BosonicObservable : public Observable {
public:
    BosonicObservable(Params& param_obj, int _freq, const std::string& _out_unit, int this_bead, BosonicExchangeBase& bosonic_exchange);

    void calculate() override;
private:
    bool bosonic;
    BosonicExchangeBase& bosonic_exchange;
};