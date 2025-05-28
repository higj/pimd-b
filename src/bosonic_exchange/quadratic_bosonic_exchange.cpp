#include <array>
#include <fstream>
#include <cmath>

#include "bosonic_exchange/quadratic_bosonic_exchange.h"

BosonicExchange::BosonicExchange(const BosonicExchangeContext& context) : BosonicExchangeBase(context),
                                                                          E_kn(context.nbosons * (context.nbosons + 1) /
                                                                              2),
                                                                          m_prefix_pot(context.nbosons + 1),
                                                                          m_suffix_pot(context.nbosons + 1),
                                                                          m_connection_probabilities(
                                                                              context.nbosons * (context.nbosons)),
                                                                          prim_est(context.nbosons + 1),
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

/*void BosonicExchange::evaluateCycleEnergies()
{
    const double half_spring_k = 0.5 * m_context.spring_constant;

    dVec x_first_bead(m_context.nbosons);
    dVec x_last_bead(m_context.nbosons);

    assignFirstLast(x_first_bead, x_last_bead);

    // Corresponds to line 5 in Algorithm 1 in the paper
    std::vector<double> diagonal_energies(m_context.nbosons);
    for (int i = 0; i < m_context.nbosons; i++)
    {
        // m_temp_nbosons_array[i] is E^[i,i]
        m_temp_nbosons_array[i] = getBeadsSeparationSquared(x_first_bead, i, x_last_bead, i);
    }

    for (int v = 0; v < m_context.nbosons; v++)
    {
        setEnk(v + 1, 1, half_spring_k * m_temp_nbosons_array[v]);

        // Corresponds to lines 6-7 in Algorithm 1 in the paper
        for (int u = v - 1; u >= 0; u--)
        {
            const double cycle_energy_of_u_to_v = getEnk(v + 1, v - u) +
                half_spring_k * (
                    // connect u to u+1
                    + getBeadsSeparationSquared(x_last_bead, u, x_first_bead, u + 1)
                    // break cycle [u+1,v]
                    - getBeadsSeparationSquared(x_first_bead, u + 1, x_last_bead, v)
                    // close cycle from v to u
                    + getBeadsSeparationSquared(x_first_bead, u, x_last_bead, v)
                    );

            setEnk(v + 1, v - u + 1, cycle_energy_of_u_to_v);
        }
    }
}*/

void BosonicExchange::evaluateCycleEnergies() {
    const double half_spring_k = 0.5 * m_context.spring_constant;

    // Get first and last bead positions
    //dVec x_first_bead(m_context.nbosons);
    //dVec x_last_bead(m_context.nbosons);
    //assignFirstLast(x_first_bead, x_last_bead);

    // Compute cycle energies using Eqs. 5-7 from the paper
    for (int v = 0; v < m_context.nbosons; ++v) {
        // Initialize single-particle cycle energy E^[v,v] (Algorithm 1, line 5)
        //const double diagonal_energy = getBeadsSeparationSquared(*m_context.x_first_bead, v, *m_context.x_last_bead, v);
        const double diagonal_energy = getExteriorSeparationSquared(v, v);
        setEnk(v + 1, 1, half_spring_k * diagonal_energy);

        // Build up multi-particle cycles (Algorithm 1, lines 6-7)
        for (int u = v - 1; u >= 0; --u) {
            const int cycle_length = v - u + 1;

            // Calculate energies needed to compute the cycle E^[u,v] from E^[u+1,v]
            //const double connect_energy = getBeadsSeparationSquared(*m_context.x_last_bead, u, *m_context.x_first_bead, u + 1);
            //const double break_energy = getBeadsSeparationSquared(*m_context.x_first_bead, u + 1, *m_context.x_last_bead, v);
            //const double close_energy = getBeadsSeparationSquared(*m_context.x_first_bead, u, *m_context.x_last_bead, v);
            const double connect_energy = getExteriorSeparationSquared(u, u + 1);
            const double break_energy = getExteriorSeparationSquared(u + 1, v);
            const double close_energy = getExteriorSeparationSquared(u, v);

            // Compute E^[u,v]
            const double previous_cycle_energy = getEnk(v + 1, v - u);
            const double spring_correction = half_spring_k * (
                connect_energy  // connect u to u+1
                - break_energy  // break cycle [u+1,v] 
                + close_energy  // close cycle from v to u
                );

            const double total_cycle_energy = previous_cycle_energy + spring_correction;

            // Validate energy before storing
            if (!std::isfinite(total_cycle_energy)) {
                throw std::runtime_error(
                    std::format("Non-finite cycle energy computed for cycle [{},{}]: previous={:e}, correction={:e}",
                        u, v, previous_cycle_energy, spring_correction)
                );
            }

            setEnk(v + 1, cycle_length, total_cycle_energy);
        }
    }
}

double BosonicExchange::getEnk(int m, int k) const
{
    const int end_of_m = m * (m + 1) / 2;
    return E_kn[end_of_m - k];
}

void BosonicExchange::setEnk(int m, int k, double val)
{
    const int end_of_m = m * (m + 1) / 2;
    E_kn[end_of_m - k] = val;
}

/*void BosonicExchange::evaluatePrefixPotential()
{
    m_prefix_pot[0] = 0.0;

    // Corresponds to lines 10-11 in Algorithm 1 in the paper
    for (int m = 1; m < m_context.nbosons + 1; m++)
    {
        double e_shift = std::numeric_limits<double>::max();

        // Shift the energies in the exponents to avoid numerical instability (Xiong & Xiong method)
        for (int k = m; k > 0; k--)
        {
            double val = getEnk(m, k) + m_prefix_pot[m - k];
            e_shift = std::min(e_shift, val);
            m_temp_nbosons_array[k - 1] = val;
        }

        double sig_denom = 0.0;
        for (int k = m; k > 0; k--)
        {
            sig_denom += exp(-m_context.thermo_beta * (m_temp_nbosons_array[k - 1] - e_shift));
        }

        m_prefix_pot[m] = e_shift - (1.0 / m_context.thermo_beta) * log(sig_denom / static_cast<double>(m));

        if (!std::isfinite(m_prefix_pot[m]))
        {
            throw std::overflow_error(
                std::format("Invalid sig_denom {:4.2f} with e_shift {:4.2f} in bosonic exchange potential", sig_denom,
                            e_shift)
            );
        }
    }
}*/

void BosonicExchange::evaluatePrefixPotential() {
    m_prefix_pot[0] = 0.0;
    std::vector<double> subdivision_potentials(m_context.nbosons);

    // Corresponds to lines 10-11 in Algorithm 1 in the paper ("last_idx" is "v" in the paper)
    for (int last_idx = 1; last_idx <= m_context.nbosons; last_idx++) {
        double e_shift = std::numeric_limits<double>::max();

        // Calculate the cycle energy of the last k particles in the sequence 1,...,last_idx
        // and add the prefix potential of the remaining particles.
        for (int k = 1; k <= last_idx; k++) {
            subdivision_potentials[k - 1] = getEnk(last_idx, k) + m_prefix_pot[last_idx - k];

            // Shift for the energies in the exponents to avoid numerical instability (Xiong & Xiong method)
            e_shift = std::min(e_shift, subdivision_potentials[k - 1]);
        }

        // Calculate the sum of the exponentials
        double sig_denom = 0.0;
        for (int k = 1; k <= last_idx; k++) {
            sig_denom += std::exp(-m_context.thermo_beta * (subdivision_potentials[k - 1] - e_shift));
        }

        // Calculate the prefix potential
        m_prefix_pot[last_idx] = e_shift - std::log(sig_denom / static_cast<double>(last_idx)) / m_context.thermo_beta;

        if (!std::isfinite(m_prefix_pot[last_idx])) {
            throw std::overflow_error(
                std::format("Invalid prefix potential calculation at last_idx={}: sig_denom={:.6e}, e_shift={:.6e}, result={:.6e}",
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

        if (!std::isfinite(m_suffix_pot[first_idx])) {
            throw std::overflow_error(
                std::format("Invalid suffix potential calculation at first_idx={}: sig_denom={:.6e}, e_shift={:.6e}, result={:.6e}",
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
    return E_kn[i];
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
            //double diff_next[NDIM];
            //getBeadsSeparation(*m_context.x_last_bead, l, *m_context.x_first_bead, next_l, diff_next);

            std::array<double, NDIM> diff_next;
            getExteriorBeadsSeparation(next_l, l, diff_next);

            const double prob = m_connection_probabilities[m_context.nbosons * l + next_l];

            for (int axis = 0; axis < NDIM; ++axis)
            {
                sums[axis] += prob * diff_next[axis];
            }
        }

        /*
        double diff_prev[NDIM];
        getBeadsSeparation(*m_context.x_last_bead, l, *m_context.x_neighbor_bead, l, diff_prev);

        for (int axis = 0; axis < NDIM; ++axis)
        {
            sums[axis] += diff_prev[axis];
        }
        */

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
            //double diff_prev[NDIM];
            //getBeadsSeparation(*m_context.x_first_bead, l, *m_context.x_last_bead, prev_l, diff_prev);

            std::array<double, NDIM> diff_prev;
            getExteriorBeadsSeparation(l, prev_l, diff_prev);

            const double prob = m_connection_probabilities[m_context.nbosons * prev_l + l];

            for (int axis = 0; axis < NDIM; ++axis)
            {
                //sums[axis] += prob * diff_prev[axis];
                sums[axis] -= prob * diff_prev[axis];
            }
        }

        /*
        double diff_next[NDIM];
        getBeadsSeparation(*m_context.x_first_bead, l, *m_context.x_neighbor_bead, l, diff_next);

        for (int axis = 0; axis < NDIM; ++axis)
        {
            sums[axis] += diff_next[axis];
        }
        */

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
    return exp(-m_context.thermo_beta * (getEnk(m_context.nbosons, m_context.nbosons) - m_prefix_pot[m_context.nbosons]));
}

double BosonicExchange::primitiveEnergyEstimator()
{
    /// TODO: Perhaps move the definition of prim_est to this method instead of having it as member variable?
    prim_est[0] = 0.0;

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

            sig += (prim_est[m - k] - e_kn_val) * exp(-m_context.thermo_beta * (e_kn_val + m_prefix_pot[m - k] - e_shift));
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
                debug << std::format("P[l={}, u={}] = {}\n", l, u, m_connection_probabilities[m_context.nbosons * l + u]);
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
