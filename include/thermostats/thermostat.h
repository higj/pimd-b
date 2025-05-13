#pragma once

#include <vector>
#include <memory>
#include "params.h"

class Coupling;

class Thermostat {
public:
    explicit Thermostat(Coupling& coupling, Params& param_obj);
    virtual ~Thermostat() = default;
    void step();
    virtual void momentaUpdate();
    virtual double getAdditionToH();
protected:
    Coupling& coupling;
    int natoms;
    double mass;
    double dt;
    int nbeads;
    double beta;
    double thermo_beta;
};