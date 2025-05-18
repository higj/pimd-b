#pragma once

#include <memory>

class RandomGenerators;

/**
 * Parameters unique to the Langevin thermostat.
 */
struct LangevinThermostatContext {
    // Generator for the random noise
    std::shared_ptr<RandomGenerators> rng;

    // Friction coefficient
    double gamma;
};