#pragma once

#include "observables/observable.h"

class Simulation; // Forward declaration

/* -------------------------------- */

class BosonicObservable : public Observable {
public:
    BosonicObservable(const Simulation& _sim, int _freq, const std::string& _out_unit);

    void calculate() override;
};