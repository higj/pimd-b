#include "states.h"
#include "simulation.h"

/**
 * @brief Generic state class constructor
 */
State::State(Params& param_obj, int _freq, const std::string& _out_unit) :
    freq(_freq), out_unit(_out_unit) {
    getVariant(param_obj.sim["nbeads"], nbeads);
    getVariant(param_obj.sys["natoms"], natoms);
}

/**
 * @brief Close the file upon destruction.
 */
State::~State() {
    if (out_file.is_open()) {
        out_file.close();
    }
}

/**
 * @brief Creates a state object based on the state type.
 *
 * @param state_type Type of state to create
 * @param _sim Simulation object
 * @param _freq Frequency at which the state is recorded
 * @param _out_unit Output unit of the state
 * @return std::unique_ptr<State> Pointer to the created state object
 * @throw std::invalid_argument If the state type is unknown
 */
std::unique_ptr<State> StateFactory::createQuantity(const std::string& state_type,
    Simulation& sim, 
    Params& param_obj,
    int _freq,
    const std::string& _out_unit) {
    if (state_type == "position") {
        return std::make_unique<PositionState>(param_obj, _freq, _out_unit, sim.coord);
    } else if (state_type == "velocity") {
        return std::make_unique<VelocityState>(param_obj, _freq, _out_unit, sim.momenta);
    } else if (state_type == "force") {
        return std::make_unique<ForceState>(param_obj, _freq, _out_unit, sim.forces);
    } else {
        throw std::invalid_argument("Unknown state type.");
    }
}