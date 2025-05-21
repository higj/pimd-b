#pragma once

#include "observables/energy.h"

class Potential;

/* -------------------------------- */

class QuantumObservable : public EnergyObservable {
public:
    QuantumObservable(Params& param_obj, int _freq, const std::string& _out_unit, int this_bead, 
                      dVec& prev_coord, dVec& coord, BosonicExchangeBase& bosonic_exchange,
                      Potential& ext_potential, Potential& int_potential, dVec& physical_forces, int& md_step);

    void calculate() override;

private:
    void calculateKinetic();
    void calculatePotential();
    std::string external_potential_name;
    std::string interaction_potential_name;
    double beta;
    Potential& ext_potential;
    Potential& int_potential;
    dVec& physical_forces;
    int& md_step;
};