#include "bosonic_exchange/fictitious_exchange.h"

#include <array>
#include <cmath>

FictitiousExchange::FictitiousExchange(
    const std::shared_ptr<const VecArray>& coord_first_bead,
    const std::shared_ptr<const VecArray>& coord_last_bead,
    const ThermalContext& thermal_ctx,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    const BeadContext& bead_ctx,
    double exchange_xi
) : BosonicExchangeBase(
    coord_first_bead,
    coord_last_bead,
    thermal_ctx,
    spring_ctx,
    box_ctx,
    bead_ctx),
    m_xi(exchange_xi),
    m_cycle_energies(bead_ctx.natoms* (bead_ctx.natoms + 1) / 2),
    m_prefix_log_absw(bead_ctx.natoms + 1),
    m_prefix_sign(bead_ctx.natoms + 1),
    m_suffix_log_absw(bead_ctx.natoms + 1),
    m_suffix_sign(bead_ctx.natoms + 1),
    m_connection_factors(bead_ctx.natoms* (bead_ctx.natoms)),
    m_log_n_factorial(std::lgamma(bead_ctx.natoms + 1))
{
    evaluateBosonicEnergies();
}

void FictitiousExchange::evaluateBosonicEnergies() {
    evaluateCycleEnergies();
    evaluatePrefixWeight();
    evaluateSuffixWeight();
    evaluateConnectionFactors();
}

void FictitiousExchange::prepare() {
    evaluateBosonicEnergies();
}

void FictitiousExchange::evaluateCycleEnergies() {
    const double half_spring_k = 0.5 * m_spring_ctx.spring_constant;

    // Compute cycle energies using Eqs. 5-7 from the paper
    for (int v = 0; v < m_bead_ctx.natoms; ++v) {
        // Initialize single-particle cycle energy E^[v,v] (Algorithm 1, line 5)
        const double diagonal_energy = getExteriorSeparationSquared(v, v);

        // By definition, E(k,m) = E^[m-k+1,m], so E^[v,v] corresponds to E(k=1,m=v).
        // Since we are using a 0-based index, we need to set E^[v+1,v+1] = E(1,v+1).
        setEnk(v + 1, 1, half_spring_k * diagonal_energy);

        // Build up multi-particle cycles (Algorithm 1, lines 6-7)
        for (int u = v - 1; u >= 0; --u) {
            const int cycle_length = v - u + 1;

            // Calculate squared distances needed to compute the cycle energy E^[u,v] from E^[u+1,v]
            const double connect_diff = getExteriorSeparationSquared(u + 1, u);
            const double break_diff = getExteriorSeparationSquared(u + 1, v);
            const double close_diff = getExteriorSeparationSquared(u, v);

            // Compute E^[u,v]
            const double previous_cycle_energy = getEnk(v + 1, v - u);
            // E^[u+1,v] corresponds to E(k=v-u,m=v+1), for 0-based indexing
            const double spring_correction = half_spring_k * (
                connect_diff // connect u to u+1
                - break_diff // break cycle [u+1,v] 
                + close_diff // close cycle from v to u
                );

            const double total_cycle_energy = previous_cycle_energy + spring_correction;

            // Validate energy before storing
            if (!std::isfinite(total_cycle_energy)) {
                throw std::runtime_error(
                    std::format("Non-finite cycle energy computed for cycle [{},{}]: previous={:e}, correction={:e}",
                        u, v, previous_cycle_energy, spring_correction)
                );
            }

            // Set the energies of all the cycles for the current v
            setEnk(v + 1, cycle_length, total_cycle_energy);
        }
    }
}

double FictitiousExchange::getEnk(int m, int k) const {
    const int end_of_m = m * (m + 1) / 2;
    return m_cycle_energies[end_of_m - k];
}

void FictitiousExchange::setEnk(int m, int k, double val) {
    const int end_of_m = m * (m + 1) / 2;
    m_cycle_energies[end_of_m - k] = val;
}

void FictitiousExchange::evaluatePrefixWeight() {
    const double log_abs_xi = std::log(std::abs(m_xi));
    const double sign_xi = (m_xi > 0.0) ? 1.0 : -1.0;

    m_prefix_log_absw[0] = 0.0;
    m_prefix_sign[0] = 1.0;

    std::vector<double> log_mag(m_bead_ctx.natoms);
    std::vector<double> term_sign(m_bead_ctx.natoms);

    for (int last_idx = 1; last_idx <= m_bead_ctx.natoms; last_idx++) {
        double shift = -std::numeric_limits<double>::infinity();

        for (int k = 1; k <= last_idx; k++) {
            const double prev_log_absw = m_prefix_log_absw[last_idx - k];
            const double prev_sign = m_prefix_sign[last_idx - k];

            // xi^(k-1) contributes (k-1)*log|xi| in magnitude and sign_xi^(k-1) in sign
            const double sign_xi_pow = (k % 2 == 1) ? 1.0 : sign_xi; // sign_xi^(k-1)

            log_mag[k - 1] = (k - 1) * log_abs_xi + prev_log_absw
                - m_thermal_ctx.thermo_beta * getEnk(last_idx, k);
            term_sign[k - 1] = sign_xi_pow * prev_sign;

            if (term_sign[k - 1] != 0.0)
                shift = std::max(shift, log_mag[k - 1]);
        }

        // All contributing terms vanished (rare, but possible): result is exactly zero
        if (!std::isfinite(shift)) {
            m_prefix_log_absw[last_idx] = -std::numeric_limits<double>::infinity();
            m_prefix_sign[last_idx] = 0.0;
            continue;
        }

        double signed_sum = 0.0;
        for (int k = 1; k <= last_idx; k++) {
            if (term_sign[k - 1] == 0.0) continue;
            signed_sum += term_sign[k - 1] * std::exp(log_mag[k - 1] - shift);
        }

        // Divide by N: subtract log(N) in log-space, sign unaffected
        if (signed_sum == 0.0) {
            m_prefix_log_absw[last_idx] = -std::numeric_limits<double>::infinity();
            m_prefix_sign[last_idx] = 0.0;
        } else {
            m_prefix_sign[last_idx] = (signed_sum > 0.0) ? 1.0 : -1.0;
            m_prefix_log_absw[last_idx] = shift + std::log(std::abs(signed_sum))
                - std::log(static_cast<double>(last_idx));
        }

        if (m_prefix_sign[last_idx] != 0.0 && !std::isfinite(m_prefix_log_absw[last_idx])) {
            throw std::overflow_error(
                std::format(
                    "Invalid prefix weight calculation at last_idx={}: signed_sum={:.6e}, shift={:.6e}",
                    last_idx, signed_sum, shift)
            );
        }
    }
}

/*
void FictitiousExchange::evaluateSuffixWeight() {
    const double log_abs_xi = std::log(std::abs(m_xi));
    const double sign_xi = (m_xi > 0.0) ? 1.0 : -1.0;

    // Store suffix weight components in log/sign form
    m_suffix_log_absw.resize(m_bead_ctx.natoms + 1);
    m_suffix_sign.resize(m_bead_ctx.natoms + 1);

    m_suffix_log_absw[m_bead_ctx.natoms] = 0.0;
    m_suffix_sign[m_bead_ctx.natoms] = 1.0;

    std::vector<double> log_mag(m_bead_ctx.natoms);
    std::vector<double> term_sign(m_bead_ctx.natoms);

    // Corresponds to lines 14-15 in Algorithm 1 in the paper ("first_idx" is "u" in the paper)
    for (int first_idx = m_bead_ctx.natoms - 1; first_idx >= 0; first_idx--) {
        double shift = -std::numeric_limits<double>::infinity();

        // Calculate the cycle energy of the first "ell" particles in the sequence u,...,N
        for (int ell = first_idx; ell < m_bead_ctx.natoms; ell++) {
            const double prev_log_absw = m_suffix_log_absw[ell + 1];
            const double prev_sign = m_suffix_sign[ell + 1];

            // xi^(ell-u)/ell contributes (ell-u)*log|xi| - log(ell) in magnitude
            const double sign_xi_pow = ((ell - first_idx) % 2 == 0) ? 1.0 : sign_xi; // sign_xi^(ell-u)

            log_mag[ell] = (ell - first_idx) * log_abs_xi - std::log(static_cast<double>(ell + 1))
                + prev_log_absw
                - m_thermal_ctx.thermo_beta * getEnk(ell + 1, ell - first_idx + 1);
            term_sign[ell] = sign_xi_pow * prev_sign;

            if (term_sign[ell] != 0.0)
                shift = std::max(shift, log_mag[ell]);
        }

        // All contributing terms vanished
        if (!std::isfinite(shift)) {
            m_suffix_log_absw[first_idx] = -std::numeric_limits<double>::infinity();
            m_suffix_sign[first_idx] = 0.0;
            continue;
        }

        double signed_sum = 0.0;
        for (int ell = first_idx; ell < m_bead_ctx.natoms; ell++) {
            if (term_sign[ell] == 0.0) continue;
            signed_sum += term_sign[ell] * std::exp(log_mag[ell] - shift);
        }

        if (signed_sum == 0.0) {
            m_suffix_log_absw[first_idx] = -std::numeric_limits<double>::infinity();
            m_suffix_sign[first_idx] = 0.0;
        } else {
            m_suffix_sign[first_idx] = (signed_sum > 0.0) ? 1.0 : -1.0;
            m_suffix_log_absw[first_idx] = shift + std::log(std::abs(signed_sum));
        }

        if (m_suffix_sign[first_idx] != 0.0 && !std::isfinite(m_suffix_log_absw[first_idx])) {
            throw std::overflow_error(
                std::format(
                    "Invalid suffix weight calculation at first_idx={}: signed_sum={:.6e}, shift={:.6e}",
                    first_idx, signed_sum, shift)
            );
        }
    }

    // Consistency check: first suffix should match last prefix
    if (std::abs(m_suffix_log_absw[0] - m_prefix_log_absw[m_bead_ctx.natoms]) > 1e-10) {
        throw std::logic_error(
            std::format("Prefix-suffix mismatch: prefix={:.6e}, suffix={:.6e}",
                m_prefix_log_absw[m_bead_ctx.natoms], m_suffix_log_absw[0])
        );
    }
}
*/

void FictitiousExchange::evaluateSuffixWeight() {
    const double log_abs_xi = std::log(std::abs(m_xi));
    const double sign_xi = (m_xi > 0.0) ? 1.0 : -1.0;

    m_suffix_log_absw.resize(m_bead_ctx.natoms + 1);
    m_suffix_sign.resize(m_bead_ctx.natoms + 1);
    m_suffix_log_absw[m_bead_ctx.natoms] = 0.0;
    m_suffix_sign[m_bead_ctx.natoms] = 1.0;

    std::vector<double> log_mag(m_bead_ctx.natoms);
    std::vector<double> term_sign(m_bead_ctx.natoms);

    for (int first_idx = m_bead_ctx.natoms - 1; first_idx >= 0; first_idx--) {
        double shift = -std::numeric_limits<double>::infinity();

        for (int ell = first_idx; ell < m_bead_ctx.natoms; ell++) {
            const double prev_log_absw = m_suffix_log_absw[ell + 1];
            const double prev_sign = m_suffix_sign[ell + 1];
            const double sign_xi_pow = ((ell - first_idx) % 2 == 0) ? 1.0 : sign_xi;

            // per-term divisor 1/ell
            log_mag[ell] = (ell - first_idx) * log_abs_xi - std::log(static_cast<double>(ell + 1))
                + prev_log_absw
                - m_thermal_ctx.thermo_beta * getEnk(ell + 1, ell - first_idx + 1);
            term_sign[ell] = sign_xi_pow * prev_sign;

            if (term_sign[ell] != 0.0)
                shift = std::max(shift, log_mag[ell]);
        }

        if (!std::isfinite(shift)) {
            m_suffix_log_absw[first_idx] = -std::numeric_limits<double>::infinity();
            m_suffix_sign[first_idx] = 0.0;
            continue;
        }

        double signed_sum = 0.0;
        for (int ell = first_idx; ell < m_bead_ctx.natoms; ell++) {
            if (term_sign[ell] == 0.0) continue;
            signed_sum += term_sign[ell] * std::exp(log_mag[ell] - shift);
        }

        if (signed_sum == 0.0) {
            m_suffix_log_absw[first_idx] = -std::numeric_limits<double>::infinity();
            m_suffix_sign[first_idx] = 0.0;
        } else {
            m_suffix_sign[first_idx] = (signed_sum > 0.0) ? 1.0 : -1.0;
            m_suffix_log_absw[first_idx] = shift + std::log(std::abs(signed_sum));
            // no further division - every term already carries its own 1/ell
        }

        if (m_suffix_sign[first_idx] != 0.0 && !std::isfinite(m_suffix_log_absw[first_idx])) {
            throw std::overflow_error(
                std::format(
                    "Invalid suffix weight calculation at first_idx={}: signed_sum={:.6e}, shift={:.6e}",
                    first_idx, signed_sum, shift)
            );
        }
    }

    // Consistency check: sign AND magnitude must agree; guard the -inf - (-inf) = NaN case
    const bool both_finite = std::isfinite(m_suffix_log_absw[0])
        && std::isfinite(m_prefix_log_absw[m_bead_ctx.natoms]);
    const bool magnitude_ok = both_finite &&
        std::abs(m_suffix_log_absw[0] - m_prefix_log_absw[m_bead_ctx.natoms]) < 1e-10;
    const bool zero_case_ok = !both_finite
        && m_suffix_sign[0] == 0.0 && m_prefix_sign[m_bead_ctx.natoms] == 0.0;
    const bool sign_ok = m_suffix_sign[0] == m_prefix_sign[m_bead_ctx.natoms];

    if (!sign_ok || !(magnitude_ok || zero_case_ok)) {
        throw std::logic_error(
            std::format("Prefix-suffix mismatch: prefix=({:.6e},{}), suffix=({:.6e},{})",
                m_prefix_log_absw[m_bead_ctx.natoms], m_prefix_sign[m_bead_ctx.natoms],
                m_suffix_log_absw[0], m_suffix_sign[0])
        );
    }
}

double FictitiousExchange::effectivePotential() {
    return (-1.0 / m_thermal_ctx.thermo_beta) * m_prefix_log_absw[m_bead_ctx.natoms];
}

double FictitiousExchange::getVn(int n) const {
    return (-1.0 / m_thermal_ctx.thermo_beta) * m_prefix_log_absw[n];
}

double FictitiousExchange::getEknSerialOrder(int i) const {
    return m_cycle_energies[i];
}

void FictitiousExchange::evaluateConnectionFactors() {
    // NOTE: for xi < 0 these are no longer probabilities (can be negative or >1) -
    // kept the name/array to minimize churn
    const int N = m_bead_ctx.natoms;
    const double beta = m_thermal_ctx.thermo_beta;
    const double sign_xi = (m_xi > 0.0) ? 1.0 : -1.0;
    const double log_abs_xi = std::log(std::abs(m_xi));

    const double total_log_absw = m_prefix_log_absw[N];
    const double total_sign = m_prefix_sign[N];

    if (total_sign == 0.0) {
        throw std::overflow_error(
            "evaluateConnectionFactors: total prefix weight is exactly zero, "
            "connection factors are ill-defined (complete sign cancellation).");
    }

    // Corresponds to lines 16-17 in Algorithm 1 in the paper - direct continuation, no xi factor
    // of its own (it's already folded into m_prefix_log_absw[l+1] via the recursion).
    for (int l = 0; l < N - 1; l++) {
        const double logmag = m_prefix_log_absw[l + 1] + m_suffix_log_absw[l + 1] - total_log_absw;
        const double sign = m_prefix_sign[l + 1] * m_suffix_sign[l + 1] * total_sign;

        m_connection_factors[N * l + (l + 1)] = 1.0 - sign * std::exp(logmag);
    }

    // Corresponds to lines 18-20 in Algorithm 1 in the paper - G(sigma)(l) = u,
    // cycle length k = l-u+1, carries the xi^{l-u} = xi^{k-1} factor explicitly.
    for (int u = 0; u < N; u++) {
        for (int l = u; l < N; l++) {
            const int k = l - u + 1;
            const double sign_xi_pow = ((k - 1) % 2 == 0) ? 1.0 : sign_xi;

            const double logmag = m_prefix_log_absw[u]
                + (k - 1) * log_abs_xi
                - beta * getEnk(l + 1, k)
                + m_suffix_log_absw[l + 1]
                - total_log_absw
                - std::log(static_cast<double>(l + 1));
            const double sign = m_prefix_sign[u] * sign_xi_pow * m_suffix_sign[l + 1] * total_sign;

            m_connection_factors[N * l + u] = sign * std::exp(logmag);
        }
    }
}

void FictitiousExchange::springForceLastBead(VecArray& f) {
    for (int l = 0; l < m_bead_ctx.natoms; l++) {
        std::array<double, NDIM> sums = {};

        for (int next_l = 0; next_l <= l + 1 && next_l < m_bead_ctx.natoms; next_l++) {
            std::array<double, NDIM> diff_next;
            getExteriorBeadsSeparation(next_l, l, diff_next);

            const double prob = m_connection_factors[m_bead_ctx.natoms * l + next_l];

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] += prob * diff_next[axis];
            }
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) = sums[axis] * m_spring_ctx.spring_constant;
        }
    }
}

void FictitiousExchange::springForceFirstBead(VecArray& f) {
    for (int l = 0; l < m_bead_ctx.natoms; l++) {
        std::array<double, NDIM> sums = {};

        for (int prev_l = std::max(0, l - 1); prev_l < m_bead_ctx.natoms; prev_l++) {
            std::array<double, NDIM> diff_prev;
            getExteriorBeadsSeparation(l, prev_l, diff_prev);

            const double prob = m_connection_factors[m_bead_ctx.natoms * prev_l + l];

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] -= prob * diff_prev[axis];
            }
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) = sums[axis] * m_spring_ctx.spring_constant;
        }
    }
}

double FictitiousExchange::getDistinctProbability() {
    double cycle_energy_sum = 0.0;
    for (int m = 1; m <= m_bead_ctx.natoms; ++m) {
        cycle_energy_sum += getEnk(m, 1);
    }

    const int N = m_bead_ctx.natoms;
    if (m_prefix_sign[N] == 0.0) {
        throw std::overflow_error("getDistinctProbability: W(N) vanished exactly.");
    }

    const double log_val = -m_thermal_ctx.thermo_beta * cycle_energy_sum
        - m_prefix_log_absw[N] - m_log_n_factorial;
    return m_prefix_sign[N] * std::exp(log_val);
}

double FictitiousExchange::getLongestProbability() {
    const int N = m_bead_ctx.natoms;
    if (m_prefix_sign[N] == 0.0) {
        throw std::overflow_error("getLongestProbability: W(N) vanished exactly.");
    }
    const double log_val = -m_thermal_ctx.thermo_beta * getEnk(N, N) - m_prefix_log_absw[N];
    return m_prefix_sign[N] * std::exp(log_val);
}

double FictitiousExchange::primitiveEnergyEstimator() {
    const double beta = m_thermal_ctx.thermo_beta;
    const double sign_xi = (m_xi > 0.0) ? 1.0 : -1.0;
    const double log_abs_xi = std::log(std::abs(m_xi));

    std::vector<double> prim_est(m_bead_ctx.natoms + 1, 0.0); // prim_est[0] = 0 (base case)
    std::vector<double> log_mag(m_bead_ctx.natoms);
    std::vector<double> raw_sign(m_bead_ctx.natoms);
    std::vector<double> coeff(m_bead_ctx.natoms);

    for (int m = 1; m <= m_bead_ctx.natoms; ++m) {
        double shift = -std::numeric_limits<double>::infinity();

        for (int k = m; k > 0; --k) {
            const double e_kn_val = getEnk(m, k);
            const int idx = m - k;
            const double sign_xi_pow = ((k - 1) % 2 == 0) ? 1.0 : sign_xi;

            log_mag[idx] = (k - 1) * log_abs_xi + m_prefix_log_absw[m - k] - beta * e_kn_val;
            raw_sign[idx] = sign_xi_pow * m_prefix_sign[m - k];
            coeff[idx] = prim_est[m - k] - e_kn_val;

            if (raw_sign[idx] != 0.0)
                shift = std::max(shift, log_mag[idx]);
        }

        if (!std::isfinite(shift)) {
            throw std::overflow_error(
                std::format("primitiveEnergyEstimator: all terms vanished at m={}", m));
        }

        double sig = 0.0;
        for (int idx = 0; idx < m; ++idx) {
            if (raw_sign[idx] == 0.0) continue;
            sig += coeff[idx] * raw_sign[idx] * std::exp(log_mag[idx] - shift);
        }

        if (m_prefix_sign[m] == 0.0) {
            throw std::overflow_error(
                std::format("primitiveEnergyEstimator: W(m) vanished exactly at m={}", m));
        }
        const double denom = static_cast<double>(m) * m_prefix_sign[m]
            * std::exp(m_prefix_log_absw[m] - shift);

        prim_est[m] = sig / denom;

        if (!std::isfinite(prim_est[m])) {
            throw std::overflow_error(
                std::format("primitiveEnergyEstimator: non-finite result at m={}", m));
        }
    }

#if TAU_CONVENTION
    return prim_est[m_bead_ctx.natoms] / m_bead_ctx.nbeads;
#else
    return prim_est[m_bead_ctx.natoms];
#endif
}

double FictitiousExchange::getSign() {
    throw std::runtime_error("getSign() is not implemented for FictitiousExchange.");
}

void FictitiousExchange::printBosonicDebug() {

}
