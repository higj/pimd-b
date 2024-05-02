#include "bosonic_exchange.h"

BosonicExchange::BosonicExchange(int nbosons_, int np_, int bead_num_, double beta_, double spring_constant_,
                                 const dVec& x_, const dVec& x_prev_, const dVec& x_next_, bool pbc_, double size_)
    : BosonicExchangeBase(nbosons_, np_, bead_num_, beta_, spring_constant_, x_, x_prev_, x_next_, pbc_, size_),
      E_kn(nbosons_ * (nbosons_ + 1) / 2),
      V(nbosons_ + 1),
      V_backwards(nbosons_ + 1),
      connection_probabilities(nbosons_ * nbosons_),
      temp_nbosons_array(nbosons_),
      prim_est(nbosons_ + 1) {
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
    dVec x_first_bead(nbosons);
    dVec x_last_bead(nbosons);

    assignFirstLast(x_first_bead, x_last_bead);

    for (int i = 0; i < nbosons; i++) {
        // temp_nbosons_array[i] is E^[i,i]
        temp_nbosons_array[i] = getBeadsSeparationSquared(x_first_bead, i, x_last_bead, i);
    }

    for (int v = 0; v < nbosons; v++) {
        setEnk(v + 1, 1, 0.5 * spring_constant * temp_nbosons_array[v]);

        for (int u = v - 1; u >= 0; u--) {
            double val = getEnk(v + 1, v - u) +
                0.5 * spring_constant * (
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

    for (int m = 1; m < nbosons + 1; m++) {
        double e_shift = std::numeric_limits<double>::max();

        for (int k = m; k > 0; k--) {
            double val = getEnk(m, k) + V[m - k];
            e_shift = std::min(e_shift, val);
            temp_nbosons_array[k - 1] = val;
        }

        double sig_denom = 0.0;
        for (int k = m; k > 0; k--) {
            sig_denom += exp(-beta * (temp_nbosons_array[k - 1] - e_shift));
        }

        V[m] = e_shift - (1.0 / beta) * log(sig_denom / static_cast<double>(m));

        if (!std::isfinite(V[m])) {
            throw std::overflow_error(
                std::format("Invalid sig_denom {:4.2f} with e_shift {:4.2f} in bosonic exchange potential", sig_denom,
                            e_shift)
            );
        }
    }
}

void BosonicExchange::evaluateVBackwards() {
    V_backwards[nbosons] = 0.0;

    for (int l = nbosons - 1; l > 0; l--) {
        double e_shift = std::numeric_limits<double>::max();
        for (int p = l; p < nbosons; p++) {
            double val = getEnk(p + 1, p - l + 1) + V_backwards[p + 1];
            e_shift = std::min(e_shift, val);
            temp_nbosons_array[p] = val;
        }

        double sig_denom = 0.0;
        for (int p = l; p < nbosons; p++) {
            sig_denom += 1.0 / (p + 1) * exp(-beta * (temp_nbosons_array[p] - e_shift));
        }

        V_backwards[l] = e_shift - log(sig_denom) / beta;

        if (!std::isfinite(V_backwards[l])) {
            throw std::overflow_error(
                std::format("Invalid sig_denom {:4.2f} with e_shift {:4.2f} in bosonic exchange potential", sig_denom,
                            e_shift)
            );
        }
    }

    V_backwards[0] = V[nbosons];
}

double BosonicExchange::effectivePotential() {
    return V[nbosons];
}

double BosonicExchange::getVn(int n) const {
    return V[n];
}

double BosonicExchange::getEknSerialOrder(int i) const {
    return E_kn[i];
}

void BosonicExchange::evaluateConnectionProbabilities() {
    for (int l = 0; l < nbosons - 1; l++) {
        double direct_link_probability = 1.0 - (exp(-beta *
            (V[l + 1] + V_backwards[l + 1] -
                V[nbosons])));
        connection_probabilities[nbosons * l + (l + 1)] = direct_link_probability;
    }
    for (int u = 0; u < nbosons; u++) {
        for (int l = u; l < nbosons; l++) {
            double close_cycle_probability = 1.0 / (l + 1) *
                exp(-beta * (V[u] + getEnk(l + 1, l - u + 1) + V_backwards[l + 1]
                    - V[nbosons]));
            connection_probabilities[nbosons * l + u] = close_cycle_probability;
        }
    }
}

void BosonicExchange::springForceLastBead(dVec& f) {
    for (int l = 0; l < nbosons; l++) {
        std::vector<double> sums(NDIM, 0.0);

        for (int next_l = 0; next_l <= l + 1 && next_l < nbosons; next_l++) {
            double diff_next[NDIM];

            getBeadsSeparation(x, l, x_next, next_l, diff_next);

            double prob = connection_probabilities[nbosons * l + next_l];

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] += prob * diff_next[axis];
            }
        }

        double diff_prev[NDIM];
        getBeadsSeparation(x, l, x_prev, l, diff_prev);

        for (int axis = 0; axis < NDIM; ++axis) {
            sums[axis] += diff_prev[axis];
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) = sums[axis] * spring_constant;
        }
    }
}

void BosonicExchange::springForceFirstBead(dVec& f) {
    for (int l = 0; l < nbosons; l++) {
        std::vector<double> sums(NDIM, 0.0);

        for (int prev_l = std::max(0, l - 1); prev_l < nbosons; prev_l++) {
            double diff_prev[NDIM];

            getBeadsSeparation(x, l, x_prev, prev_l, diff_prev);

            double prob = connection_probabilities[nbosons * prev_l + l];

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] += prob * diff_prev[axis];
            }
        }

        double diff_next[NDIM];
        getBeadsSeparation(x, l, x_next, l, diff_next);

        for (int axis = 0; axis < NDIM; ++axis) {
            sums[axis] += diff_next[axis];
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) = sums[axis] * spring_constant;
        }
    }
}

/**
 * Primitive kinetic energy estimator for bosons.
 * Corresponds to Eqns. (4)-(5) in SI of pnas.1913365116, 
 * excluding the constant factor of d*N*P/(2*beta).
 * 
 * @return The exterior spring contribution to the overall kinetic energy.
 */
double BosonicExchange::primEstimator() {
    prim_est[0] = 0.0;

    for (int m = 1; m < nbosons + 1; ++m) {
        double sig = 0.0;

        // Shift the energies in the exponents to avoid numerical instability (Xiong & Xiong method)
        double e_shift = std::numeric_limits<double>::max();

        for (int k = m; k > 0; k--) {
            e_shift = std::min(e_shift, getEnk(m, k) + V[m - k]);
        }

        for (int k = m; k > 0; --k) {
            const double e_kn_val = getEnk(m, k);

            sig += (prim_est[m - k] - e_kn_val) * exp(-beta * (e_kn_val + V[m - k] - e_shift));
        }

        const double sig_denom_m = m * exp(-beta * (V[m] - e_shift));

        prim_est[m] = sig / sig_denom_m;
    }

#if IPI_CONVENTION
    return prim_est[nbosons] / nbeads;
#else
    return prim_est[nbosons];
#endif
}
