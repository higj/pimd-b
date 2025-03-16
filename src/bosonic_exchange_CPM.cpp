#include <array>
#include <fstream>
#include <cmath>
#include "mpi.h"
#include "bosonic_exchange_CPM.h"
#include "simulation.h"

BosonicExchangeCPM::BosonicExchangeCPM(const Simulation& _sim, const int physical_exchange_time) : 
    BosonicExchange(_sim), physical_exchange_time(physical_exchange_time) {}

void BosonicExchangeCPM::evaluateConnectionProbabilities() {
    if (sim.getStep()<physical_exchange_time)
        {
        double prob_threshold = 0.0;
        for (int l = 0; l < nbosons - 1; l++) {
            // Calculate the probability of connecting to oneself
            double close_cycle_probability = 1.0 / (l + 1) *
                    exp(-beta * (V[l] + getEnk(l + 1, 1) + V_backwards[l + 1] - V[nbosons]));
            // if (close_cycle_probability < prob_threshold) {
            //     connection_probabilities[nbosons * l + l] = close_cycle_probability;
            //     temp_nbosons_array[l] = 0;
            // }
            // else {
            connection_probabilities[nbosons * l + l] = prob_threshold;
            temp_nbosons_array[l] = (close_cycle_probability - prob_threshold) / (l + 1);
            // }
            // std::cout << std::format("l: {}, self_link_probability: {}\n", l, connection_probabilities[nbosons * l + l]);

            double direct_link_probability = 1.0 - (exp(-beta *
                (V[l + 1] + V_backwards[l + 1] - V[nbosons]))) + temp_nbosons_array[l];
            connection_probabilities[nbosons * l + (l + 1)] = direct_link_probability;
            // std::cout << std::format("l: {}, direct_link_probability: {}\n", l, direct_link_probability);
        }
        int l = nbosons - 1;
        double close_cycle_probability = 1.0 / (l + 1) *
                exp(-beta * (V[l] + getEnk(l + 1, 1) + V_backwards[l + 1] - V[nbosons]));
        // if (close_cycle_probability < prob_threshold) {
        //     connection_probabilities[nbosons * l + l] = close_cycle_probability;
        //     temp_nbosons_array[l] = 0;
        // }
        // else {
        connection_probabilities[nbosons * l + l] = prob_threshold;
        temp_nbosons_array[l] = (close_cycle_probability - prob_threshold) / l;
        // }
        // std::cout << std::format("l: {}, self_link_probability: {}\n", l, connection_probabilities[nbosons * l + l]);
        for (int u = 0; u < nbosons; u++) {
            for (int l = u + 1; l < nbosons; l++) {
                double close_cycle_probability = 1.0 / (l + 1) *
                    exp(-beta * (V[u] + getEnk(l + 1, l - u + 1) + V_backwards[l + 1]
                        - V[nbosons])) + temp_nbosons_array[l];
                connection_probabilities[nbosons * l + u] = close_cycle_probability;
                // std::cout << std::format("l: {}, u: {}, close_cycle_probability: {}\n", l, u, close_cycle_probability);
            }
        }
    }
    else BosonicExchange::evaluateConnectionProbabilities();
}