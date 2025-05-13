#pragma once

#include "observables/observable.h"
#include "common.h"
class Potential;

/* -------------------------------- */

class GSFActionObservable : public Observable {
public:
    GSFActionObservable(Params& param_obj, int _freq, const std::string& _out_unit, int this_bead,
                        Potential& ext_potential, Potential& int_potential, dVec& coord);

    void calculate() override;
private:
    double beta, spring_constant;
    Potential& ext_potential;
    Potential& int_potential;
    dVec& coord;
    int natoms, nbeads;
    double int_pot_cutoff, size;
    bool pbc;
};
