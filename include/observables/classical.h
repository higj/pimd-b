#pragma once

#include "observables/observable.h"

class Simulation; // Forward declaration

/* -------------------------------- */

class ClassicalObservable : public Observable {
public:
    ClassicalObservable(const Simulation& _sim, int _freq, const std::string& _out_unit);

    void calculate() override;

private:
    void calculateKineticEnergy();
    void calculateSpringEnergy();
};