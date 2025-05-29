#include <array>
#include <fstream>
#include <cmath>

#include "bosonic_exchange/quadratic_bosonic_exchange.h"

BosonicExchange::BosonicExchange(const BosonicExchangeContext& context) : BosonicExchangeBase(context),
                                                                          m_cycle_energies(context.nbosons * (context.nbosons + 1) /
                                                                              2),
                                                                          m_prefix_pot(context.nbosons + 1),
                                                                          m_suffix_pot(context.nbosons + 1),
                                                                          m_connection_probabilities(
                                                                              context.nbosons * (context.nbosons)),
                                                                          m_log_n_factorial(
                                                                              std::lgamma(context.nbosons + 1))
{
    evaluateBosonicEnergies();
}

void BosonicExchange::evaluateBosonicEnergies()
{
    evaluateCycleEnergies();
    evaluatePrefixPotential();
    evaluateSuffixPotential();
    evaluateConnectionProbabilities();
}

void BosonicExchange::prepare()
{
    evaluateBosonicEnergies();
}

void BosonicExchange::evaluateCycleEnergies()
{
    const double half_spring_k = 0.5 * m_context.spring_constant;

    // Compute cycle energies using Eqs. 5-7 from the paper
    for (int v = 0; v < m_context.nbosons; ++v)
    {
        // Initialize single-particle cycle energy E^[v,v] (Algorithm 1, line 5)
        const double diagonal_energy = getExteriorSeparationSquared(v, v);

        // By definition, E(k,m) = E^[m-k+1,m], so E^[v,v] corresponds to E(k=1,m=v).
        // Since we are using a 0-based index, we need to set E^[v+1,v+1] = E(1,v+1).
        setEnk(v + 1, 1, half_spring_k * diagonal_energy);

        // Build up multi-particle cycles (Algorithm 1, lines 6-7)
        for (int u = v - 1; u >= 0; --u)
        {
            const int cycle_length = v - u + 1;

            // Calculate squared distances needed to compute the cycle energy E^[u,v] from E^[u+1,v]
            const double connect_diff = getExteriorSeparationSquared(u + 1, u);
            const double break_diff = getExteriorSeparationSquared(u + 1, v);
            const double close_diff = getExteriorSeparationSquared(u, v);

            // Compute E^[u,v]
            const double previous_cycle_energy = getEnk(v + 1, v - u); // E^[u+1,v] corresponds to E(k=v-u,m=v+1), for 0-based indexing
            const double spring_correction = half_spring_k * (
                connect_diff // connect u to u+1
                - break_diff // break cycle [u+1,v] 
                + close_diff // close cycle from v to u
            );

            const double total_cycle_energy = previous_cycle_energy + spring_correction;

            // Validate energy before storing
            if (!std::isfinite(total_cycle_energy))
            {
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

double BosonicExchange::getEnk(int m, int k) const
{
    const int end_of_m = m * (m + 1) / 2;
    return m_cycle_energies[end_of_m - k];
}

void BosonicExchange::setEnk(int m, int k, double val)
{
    const int end_of_m = m * (m + 1) / 2;
    m_cycle_energies[end_of_m - k] = val;
}

void BosonicExchange::evaluatePrefixPotential()
{
    m_prefix_pot[0] = 0.0;
    std::vector<double> subdivision_potentials(m_context.nbosons);

    // Corresponds to lines 10-11 in Algorithm 1 in the paper ("last_idx" is "v" in the paper)
    for (int last_idx = 1; last_idx <= m_context.nbosons; last_idx++)
    {
        double e_shift = std::numeric_limits<double>::max();

        // Calculate the cycle energy of the last k particles in the sequence 1,...,last_idx
        // and add the prefix potential of the remaining particles.
        for (int k = 1; k <= last_idx; k++)
        {
            // E^[v-k+1,v] corresponds to E(k,v). Note that here both "k" and "v" ("last_idx") are 1-based indices.
            subdivision_potentials[k - 1] = getEnk(last_idx, k) + m_prefix_pot[last_idx - k];

            // Shift for the energies in the exponents to avoid numerical instability (Xiong & Xiong method)
            e_shift = std::min(e_shift, subdivision_potentials[k - 1]);
        }

        // Calculate the sum of the exponentials
        double sig_denom = 0.0;
        for (int k = 1; k <= last_idx; k++)
        {
            sig_denom += std::exp(-m_context.thermo_beta * (subdivision_potentials[k - 1] - e_shift));
        }

        // Calculate the prefix potential
        m_prefix_pot[last_idx] = e_shift - std::log(sig_denom / static_cast<double>(last_idx)) / m_context.thermo_beta;

        if (!std::isfinite(m_prefix_pot[last_idx]))
        {
            throw std::overflow_error(
                std::format(
                    "Invalid prefix potential calculation at last_idx={}: sig_denom={:.6e}, e_shift={:.6e}, result={:.6e}",
                    last_idx, sig_denom, e_shift, m_prefix_pot[last_idx])
            );
        }
    }
}

void BosonicExchange::evaluateSuffixPotential()
{
    m_suffix_pot[m_context.nbosons] = 0.0;
    std::vector<double> subdivision_potentials(m_context.nbosons);

    // Corresponds to lines 14-15 in Algorithm 1 in the paper ("first_idx" is "u" in the paper)
    for (int first_idx = m_context.nbosons - 1; first_idx > 0; first_idx--)
    {
        double e_shift = std::numeric_limits<double>::max();

        // Calculate the cycle energy of the first "ell" particles in the sequence u,...,N
        // and add the suffix potential of the remaining particles.
        for (int ell = first_idx; ell < m_context.nbosons; ell++)
        {
            // E^[u,ell] corresponds to E(ell-u+1,ell+1), since "u" ("first_idx") and "ell" are 0-based indices.
            // This is because the last element of the suffix potential, V^[N+1,N], is at index N, and not N+1.
            subdivision_potentials[ell] = getEnk(ell + 1, ell - first_idx + 1) + m_suffix_pot[ell + 1];

            // Shift for the energies in the exponents to avoid numerical instability (Xiong & Xiong method)
            e_shift = std::min(e_shift, subdivision_potentials[ell]);
        }

        // Calculate the sum of the exponentials
        double sig_denom = 0.0;
        for (int ell = first_idx; ell < m_context.nbosons; ell++)
        {
            sig_denom += 1.0 / (ell + 1) * exp(-m_context.thermo_beta * (subdivision_potentials[ell] - e_shift));
        }

        // Calculate the suffix potential
        m_suffix_pot[first_idx] = e_shift - log(sig_denom) / m_context.thermo_beta;

        if (!std::isfinite(m_suffix_pot[first_idx]))
        {
            throw std::overflow_error(
                std::format(
                    "Invalid suffix potential calculation at first_idx={}: sig_denom={:.6e}, e_shift={:.6e}, result={:.6e}",
                    first_idx, sig_denom, e_shift, m_suffix_pot[first_idx])
            );
        }
    }

    // The first suffix potential, V^[u=1,N], coincides with the last prefix potential, V^[1,v=N]
    m_suffix_pot[0] = m_prefix_pot[m_context.nbosons];
}

double BosonicExchange::effectivePotential()
{
    return m_prefix_pot[m_context.nbosons];
}

double BosonicExchange::getVn(int n) const
{
    return m_prefix_pot[n];
}

double BosonicExchange::getEknSerialOrder(int i) const
{
    return m_cycle_energies[i];
}

void BosonicExchange::evaluateConnectionProbabilities()
{
    // Corresponds to lines 16-17 in Algorithm 1 in the paper
    for (int l = 0; l < m_context.nbosons - 1; l++)
    {
        const double direct_link_probability = 1.0 - (exp(-m_context.thermo_beta *
            (m_prefix_pot[l + 1] + m_suffix_pot[l + 1] -
                m_prefix_pot[m_context.nbosons])));
        m_connection_probabilities[m_context.nbosons * l + (l + 1)] = direct_link_probability;
    }

    // Corresponds to lines 18-20 in Algorithm 1 in the paper
    for (int u = 0; u < m_context.nbosons; u++)
    {
        for (int l = u; l < m_context.nbosons; l++)
        {
            const double close_cycle_probability = 1.0 / (l + 1) *
                exp(-m_context.thermo_beta * (m_prefix_pot[u] + getEnk(l + 1, l - u + 1) + m_suffix_pot[l + 1]
                    - m_prefix_pot[m_context.nbosons]));
            m_connection_probabilities[m_context.nbosons * l + u] = close_cycle_probability;
        }
    }
}

void BosonicExchange::springForceLastBead(dVec& f)
{
    for (int l = 0; l < m_context.nbosons; l++)
    {
        std::array<double, NDIM> sums = {};

        for (int next_l = 0; next_l <= l + 1 && next_l < m_context.nbosons; next_l++)
        {
            std::array<double, NDIM> diff_next;
            getExteriorBeadsSeparation(next_l, l, diff_next);

            const double prob = m_connection_probabilities[m_context.nbosons * l + next_l];

            for (int axis = 0; axis < NDIM; ++axis)
            {
                sums[axis] += prob * diff_next[axis];
            }
        }

        for (int axis = 0; axis < NDIM; ++axis)
        {
            f(l, axis) = sums[axis] * m_context.spring_constant;
        }
    }
}

void BosonicExchange::springForceFirstBead(dVec& f)
{
    for (int l = 0; l < m_context.nbosons; l++)
    {
        std::array<double, NDIM> sums = {};

        for (int prev_l = std::max(0, l - 1); prev_l < m_context.nbosons; prev_l++)
        {
            std::array<double, NDIM> diff_prev;
            getExteriorBeadsSeparation(l, prev_l, diff_prev);

            const double prob = m_connection_probabilities[m_context.nbosons * prev_l + l];

            for (int axis = 0; axis < NDIM; ++axis)
            {
                sums[axis] -= prob * diff_prev[axis];
            }
        }

        for (int axis = 0; axis < NDIM; ++axis)
        {
            f(l, axis) = sums[axis] * m_context.spring_constant;
        }
    }
}

double BosonicExchange::getDistinctProbability()
{
    double cycle_energy_sum = 0.0;
    for (int m = 1; m < m_context.nbosons + 1; ++m)
    {
        cycle_energy_sum += getEnk(m, 1);
    }

    return exp(-m_context.thermo_beta * (cycle_energy_sum - m_prefix_pot[m_context.nbosons]) - m_log_n_factorial);
}

double BosonicExchange::getLongestProbability()
{
    return exp(
        -m_context.thermo_beta * (getEnk(m_context.nbosons, m_context.nbosons) - m_prefix_pot[m_context.nbosons]));
}

double BosonicExchange::primitiveEnergyEstimator()
{
    std::vector<double> prim_est(m_context.nbosons + 1);

    for (int m = 1; m < m_context.nbosons + 1; ++m)
    {
        double sig = 0.0;
        // Shift the energies in the exponents to avoid numerical instability (Xiong & Xiong method)
        double e_shift = std::numeric_limits<double>::max();

        for (int k = m; k > 0; k--)
        {
            e_shift = std::min(e_shift, getEnk(m, k) + m_prefix_pot[m - k]);
        }

        for (int k = m; k > 0; --k)
        {
            const double e_kn_val = getEnk(m, k);
            sig += (prim_est[m - k] - e_kn_val) * exp(
                -m_context.thermo_beta * (e_kn_val + m_prefix_pot[m - k] - e_shift));
        }

        const double sig_denom_m = m * exp(-m_context.thermo_beta * (m_prefix_pot[m] - e_shift));

        prim_est[m] = sig / sig_denom_m;
    }

#if IPI_CONVENTION
    return prim_est[m_context.nbosons] / m_context.nbeads;
#else
    return prim_est[m_context->nbosons];
#endif
}

void BosonicExchange::printBosonicDebug()
{
    if (m_context.this_bead == 0)
    {
        std::ofstream debug;
        debug.open(std::format("{}/bosonic_debug.log", Output::FOLDER_NAME), std::ios::out | std::ios::app);

        /// TODO: Think about access to number of steps; Or maybe move this method to Simulation class?
        //debug << "Step " << sim.getStep() << '\n';

        debug << "Bosonic energies:\n";
        for (int m = 1; m < m_context.nbosons + 1; ++m)
        {
            debug << "m_prefix_pot[" << m << "] = " << m_prefix_pot[m] << '\n';
        }

        debug << "----\n";

        debug << "Connection probabilities:\n";
        for (int l = 0; l < m_context.nbosons; ++l)
        {
            for (int u = 0; u < m_context.nbosons; ++u)
            {
                debug << std::format("P[l={}, u={}] = {}\n", l, u,
                                     m_connection_probabilities[m_context.nbosons * l + u]);
            }
        }

        debug << "----\n";

        debug << "getEnk(0, 0) = " << getEnk(0, 0) << '\n';

        for (int m = 1; m < m_context.nbosons + 1; ++m)
        {
            for (int k = m; k > 0; --k)
            {
                debug << std::format("getEnk(m = {}, k = {}) = {}\n", m, k, getEnk(m, k));
            }
        }

        debug << "===========\n";

        debug.close();
    }
}
