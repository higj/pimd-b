#pragma once

#include "observables/energy.h"

class Thermostat; // Forward declaration

/* -------------------------------- */

class ClassicalObservable : public EnergyObservable {
public:
    ClassicalObservable(Params& param_obj, int _freq, const std::string& _out_unit, int this_bead, 
                        dVec& prev_coord, dVec& coord, BosonicExchangeBase& bosonic_exchange, Thermostat& thermostat, dVec& momenta);

    void calculate() override;

private:
    void calculateKineticEnergy();
    void calculateSpringEnergy();
    std::string thermostat_type;
    Thermostat& thermostat;
    dVec& momenta;
};