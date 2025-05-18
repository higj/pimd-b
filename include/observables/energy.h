#pragma once

#include "observables/observable.h"
#include "common.h"

class BosonicExchangeBase; // Forward declaration
/* -------------------------------- */

class EnergyObservable : public Observable {
public:
    EnergyObservable(Params& param_obj, int _freq, const std::string& _out_unit, int this_bead, 
                     dVec& prev_coord, dVec& coord, BosonicExchangeBase& bosonic_exchange);

    void calculate() override;

protected:
    double classicalSpringEnergy() const;
    int nbeads;
    int natoms;
    bool bosonic;
    bool pbc;
    double spring_constant;
    double size;
    double mass;
    dVec& prev_coord;
    dVec& coord;
    BosonicExchangeBase& bosonic_exchange;
};