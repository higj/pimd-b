#include "states.h"
#include "simulation.h"

/**
 * @brief Generic state class constructor
 */
State::State(const Simulation& _sim, int _freq, const std::string& _out_unit) :
    sim(_sim), freq(_freq), out_unit(_out_unit) {
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
    const Simulation& _sim,
    int _freq,
    const std::string& _out_unit) {
    if (state_type == "position") {
        return std::make_unique<PositionState>(_sim, _freq, _out_unit);
    } else if (state_type == "velocity") {
        return std::make_unique<VelocityState>(_sim, _freq, _out_unit);
    } else if (state_type == "force") {
        return std::make_unique<ForceState>(_sim, _freq, _out_unit);
    } else {
        throw std::invalid_argument("Unknown state type.");
    }
}