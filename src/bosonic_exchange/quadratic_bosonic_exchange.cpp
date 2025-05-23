#include <array>
#include <fstream>
#include <cmath>

#include "bosonic_exchange/quadratic_bosonic_exchange.h"

BosonicExchange::BosonicExchange(const BosonicExchangeContext& context) : BosonicExchangeBase(context),
    E_kn(context.nbosons * (context.nbosons + 1) / 2),
    V(context.nbosons + 1),
    V_backwards(context.nbosons + 1),
    connection_probabilities(context.nbosons * (context.nbosons)),
    temp_nbosons_array(context.nbosons),
    prim_est(context.nbosons + 1),
    log_n_factorial(std::lgamma(context.nbosons + 1)) {
    evaluateBosonicEnergies();
}

void BosonicExchange::evaluateBosonicEnergies() {
    evaluateCycleEnergies();
    evaluateVBn();
    evaluateVBackwards();
    evaluateConnectionProbabilities();
}

/**
 * @brief Re-evaluate the bosonic energies and connection probabilities.
 * Typically used after coordinate updates.
 */
void BosonicExchange::prepare() {
    evaluateBosonicEnergies();
}

void BosonicExchange::evaluateCycleEnergies() {
    dVec x_first_bead(m_context.nbosons);
    dVec x_last_bead(m_context.nbosons);

    assignFirstLast(x_first_bead, x_last_bead);

    for (int i = 0; i < m_context.nbosons; i++) {
        // temp_nbosons_array[i] is E^[i,i]
        temp_nbosons_array[i] = getBeadsSeparationSquared(x_first_bead, i, x_last_bead, i);
    }

    for (int v = 0; v < m_context.nbosons; v++) {
        setEnk(v + 1, 1, 0.5 * m_context.spring_constant * temp_nbosons_array[v]);

        for (int u = v - 1; u >= 0; u--) {
            double val = getEnk(v + 1, v - u) +
                0.5 * m_context.spring_constant * (
                    // connect u to u+1
                    + getBeadsSeparationSquared(x_last_bead, u, x_first_bead, u + 1)
                    // break cycle [u+1,v]
                    - getBeadsSeparationSquared(x_first_bead, u + 1, x_last_bead, v)
                    // close cycle from v to u
                    + getBeadsSeparationSquared(x_first_bead, u, x_last_bead, v));

            setEnk(v + 1, v - u + 1, val);
        }
    }
}

double BosonicExchange::getEnk(int m, int k) const {
    int end_of_m = m * (m + 1) / 2;
    return E_kn[end_of_m - k];
}

void BosonicExchange::setEnk(int m, int k, double val) {
    int end_of_m = m * (m + 1) / 2;
    E_kn[end_of_m - k] = val;
}

void BosonicExchange::evaluateVBn() {
    V[0] = 0.0;

    for (int m = 1; m < m_context.nbosons + 1; m++) {
        double e_shift = std::numeric_limits<double>::max();

        for (int k = m; k > 0; k--) {
            double val = getEnk(m, k) + V[m - k];
            e_shift = std::min(e_shift, val);
            temp_nbosons_array[k - 1] = val;
        }

        double sig_denom = 0.0;
        for (int k = m; k > 0; k--) {
            sig_denom += exp(-m_context.thermo_beta * (temp_nbosons_array[k - 1] - e_shift));
        }

        V[m] = e_shift - (1.0 / m_context.thermo_beta) * log(sig_denom / static_cast<double>(m));

        if (!std::isfinite(V[m])) {
            throw std::overflow_error(
                std::format("Invalid sig_denom {:4.2f} with e_shift {:4.2f} in bosonic exchange potential", sig_denom,
                            e_shift)
            );
        }
    }
}

void BosonicExchange::evaluateVBackwards() {
    V_backwards[m_context.nbosons] = 0.0;

    for (int l = m_context.nbosons - 1; l > 0; l--) {
        double e_shift = std::numeric_limits<double>::max();
        for (int p = l; p < m_context.nbosons; p++) {
            double val = getEnk(p + 1, p - l + 1) + V_backwards[p + 1];
            e_shift = std::min(e_shift, val);
            temp_nbosons_array[p] = val;
        }

        double sig_denom = 0.0;
        for (int p = l; p < m_context.nbosons; p++) {
            sig_denom += 1.0 / (p + 1) * exp(-m_context.thermo_beta * (temp_nbosons_array[p] - e_shift));
        }

        V_backwards[l] = e_shift - log(sig_denom) / m_context.thermo_beta;

        if (!std::isfinite(V_backwards[l])) {
            throw std::overflow_error(
                std::format("Invalid sig_denom {:4.2f} with e_shift {:4.2f} in bosonic exchange potential", sig_denom,
                            e_shift)
            );
        }
    }

    V_backwards[0] = V[m_context.nbosons];
}

double BosonicExchange::effectivePotential() {
    return V[m_context.nbosons];
}

double BosonicExchange::getVn(int n) const {
    return V[n];
}

double BosonicExchange::getEknSerialOrder(int i) const {
    return E_kn[i];
}

void BosonicExchange::evaluateConnectionProbabilities() {
    for (int l = 0; l < m_context.nbosons - 1; l++) {
        double direct_link_probability = 1.0 - (exp(-m_context.thermo_beta *
            (V[l + 1] + V_backwards[l + 1] -
                V[m_context.nbosons])));
        connection_probabilities[m_context.nbosons * l + (l + 1)] = direct_link_probability;
    }
    for (int u = 0; u < m_context.nbosons; u++) {
        for (int l = u; l < m_context.nbosons; l++) {
            double close_cycle_probability = 1.0 / (l + 1) *
                exp(-m_context.thermo_beta * (V[u] + getEnk(l + 1, l - u + 1) + V_backwards[l + 1]
                    - V[m_context.nbosons]));
            connection_probabilities[m_context.nbosons * l + u] = close_cycle_probability;
        }
    }
}

void BosonicExchange::springForceLastBead(dVec& f) {
    for (int l = 0; l < m_context.nbosons; l++) {
        std::array<double, NDIM> sums = {};

        for (int next_l = 0; next_l <= l + 1 && next_l < m_context.nbosons; next_l++) {
            double diff_next[NDIM];

            getBeadsSeparation(*m_context.x, l, *m_context.x_next, next_l, diff_next);

            double prob = connection_probabilities[m_context.nbosons * l + next_l];

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] += prob * diff_next[axis];
            }
        }

        double diff_prev[NDIM];
        getBeadsSeparation(*m_context.x, l, *m_context.x_prev, l, diff_prev);

        for (int axis = 0; axis < NDIM; ++axis) {
            sums[axis] += diff_prev[axis];
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) = sums[axis] * m_context.spring_constant;
        }
    }
}

void BosonicExchange::springForceFirstBead(dVec& f) {
    for (int l = 0; l < m_context.nbosons; l++) {
        std::array<double, NDIM> sums = {};

        for (int prev_l = std::max(0, l - 1); prev_l < m_context.nbosons; prev_l++) {
            double diff_prev[NDIM];

            getBeadsSeparation(*m_context.x, l, *m_context.x_prev, prev_l, diff_prev);

            double prob = connection_probabilities[m_context.nbosons * prev_l + l];

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] += prob * diff_prev[axis];
            }
        }

        double diff_next[NDIM];
        getBeadsSeparation(*m_context.x, l, *m_context.x_next, l, diff_next);

        for (int axis = 0; axis < NDIM; ++axis) {
            sums[axis] += diff_next[axis];
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) = sums[axis] * m_context.spring_constant;
        }
    }
}

/**
 * Evaluate the probability of the configuration where all the particles are separate.
 *
 * @return Probability of a configuration corresponding to the identity permutation.
 */
double BosonicExchange::getDistinctProbability() {
    double cycle_energy_sum = 0.0;
    for (int m = 1; m < m_context.nbosons + 1; ++m) {
        cycle_energy_sum += getEnk(m, 1);
    }

    return exp(-m_context.thermo_beta * (cycle_energy_sum - V[m_context.nbosons]) - log_n_factorial);
}

/**
 * Evaluate the probability of a configuration where all the particles are connected,
 * divided by 1/N. Notice that there are (N-1)! permutations of this topology
 * (all represented by the cycle 0,1,...,N-1,0); this cancels the division by 1/N.
 *
 * @return Probability of a configuration where all the particles are connected.
 */
double BosonicExchange::getLongestProbability() {
    return exp(-m_context.thermo_beta * (getEnk(m_context.nbosons, m_context.nbosons) - V[m_context.nbosons]));
}

/**
 * Evaluates the partial derivative of (beta*V_B) with respect to beta.
 * This is used to calculate the primitive kinetic energy estimator for bosons.
 * Corresponds to Eqns. (4)-(5) in SI of pnas.1913365116,
 * excluding the constant factor of d*N*P/(2*beta).
 *
 * @return The exterior spring contribution to the overall kinetic energy.
 */
double BosonicExchange::primEstimator() {
    prim_est[0] = 0.0;

    for (int m = 1; m < m_context.nbosons + 1; ++m) {
        double sig = 0.0;

        // Shift the energies in the exponents to avoid numerical instability (Xiong & Xiong method)
        double e_shift = std::numeric_limits<double>::max();

        for (int k = m; k > 0; k--) {
            e_shift = std::min(e_shift, getEnk(m, k) + V[m - k]);
        }

        for (int k = m; k > 0; --k) {
            const double e_kn_val = getEnk(m, k);

            sig += (prim_est[m - k] - e_kn_val) * exp(-m_context.thermo_beta * (e_kn_val + V[m - k] - e_shift));
        }

        const double sig_denom_m = m * exp(-m_context.thermo_beta * (V[m] - e_shift));

        prim_est[m] = sig / sig_denom_m;
    }

#if IPI_CONVENTION
    return prim_est[m_context.nbosons] / m_context.nbeads;
#else
    return prim_est[m_context->nbosons];
#endif
}

void BosonicExchange::printBosonicDebug() {
    if (m_context.this_bead == 0) {
        std::ofstream debug;
        debug.open(std::format("{}/bosonic_debug.log", Output::FOLDER_NAME), std::ios::out | std::ios::app);

        /// TODO: Think about access to number of steps; Or maybe move this method to Simulation class?
        //debug << "Step " << sim.getStep() << '\n';

        debug << "Bosonic energies:\n";
        for (int m = 1; m < m_context.nbosons + 1; ++m) {
            debug << "V[" << m << "] = " << V[m] << '\n';
        }

        debug << "----\n";

        debug << "Connection probabilities:\n";
        for (int l = 0; l < m_context.nbosons; ++l) {
            for (int u = 0; u < m_context.nbosons; ++u) {
                debug << std::format("P[l={}, u={}] = {}\n", l, u, connection_probabilities[m_context.nbosons * l + u]);
            }
        }

        debug << "----\n";

        debug << "getEnk(0, 0) = " << getEnk(0, 0) << '\n';

        for (int m = 1; m < m_context.nbosons + 1; ++m) {
            for (int k = m; k > 0; --k) {
                debug << std::format("getEnk(m = {}, k = {}) = {}\n", m, k, getEnk(m, k));
            }
        }

        debug << "===========\n";

        debug.close();
    }
}