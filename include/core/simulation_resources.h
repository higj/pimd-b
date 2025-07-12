#pragma once

#include <memory>
#include <vector>

struct SimulationConfig;
class SystemState;
class ExchangeState;
class RandomGenerators;
class ForceManager;
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
    std::shared_ptr<ForceManager> force_mgr;

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
