#pragma once

#include <memory>
#include <vector>
#include "simulation_config.h"

class SystemState;
struct ExchangeState;
class RandomGenerators;
//class Box;
class ForceFieldManager;
class Propagator;
class Thermostat;
class NormalModes;
class Observable;
class Dump;

struct SimulationResources {
    std::shared_ptr<SimulationConfig> config;
    std::shared_ptr<SystemState> state;
    std::shared_ptr<ExchangeState> exchange_state;
    std::shared_ptr<RandomGenerators> rng;
    //std::shared_ptr<Box> box;
    std::shared_ptr<ForceFieldManager> force_mgr;

    std::shared_ptr<NormalModes> normal_modes;
    std::shared_ptr<Propagator> propagator;
    std::shared_ptr<Thermostat> thermostat;
    /*struct ThermostatResources {
        std::shared_ptr<Thermostat> thermostat;
        VariantMap thermostat_params;
    };*/

    std::vector<std::shared_ptr<Observable>> observables;
    std::vector<std::shared_ptr<Dump>> dumps;
};
