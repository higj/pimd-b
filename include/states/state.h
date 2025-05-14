#pragma once

#include <string>
#include <memory>
#include <fstream>
#include "common.h"

class Params; // Forward declaration
class Simulation; // Forward declaration

/* -------------------------------- */

class State {
public:
    explicit State(Params& param_obj, int _freq, const std::string& _out_unit);
    ~State();

    virtual void initialize(int this_bead) = 0;
    virtual void output(int step) = 0;   

protected:
    int freq;                // Frequency at which the state is recorded
    std::ofstream out_file;  // Output file stream
    std::string out_unit;    // Units of the output quantities
    int natoms;
};

/* -------------------------------- */

class StateFactory {
public:
    static std::unique_ptr<State> createQuantity(const std::string& state_type, 
        Simulation& _sim, Params& param_obj, int _freq, const std::string& _out_unit);
};