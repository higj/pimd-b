#pragma once

#include "observables/observable.h"

class Simulation; // Forward declaration

/* -------------------------------- */

class EnergyObservable : public Observable {
public:
    EnergyObservable(const Simulation& _sim, int _freq, const std::string& _out_unit);

    void calculate() override;

private:
    void calculateKinetic();
    void calculatePotential();
};