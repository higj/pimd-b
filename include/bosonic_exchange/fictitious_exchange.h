#pragma once

#include "bosonic_exchange/bosonic_exchange_base.h"

#include <vector>

/**
 * @class FictitiousExchange
 * @brief Implements the quadratic-scaling bosonic exchange algorithm for path integral molecular dynamics.
 *
 * This class implements the bosonic exchange algorithm by Yotam M. Y. Feldman and Barak Hirshberg (DOI: 10.1063/5.0173749).
 * It manages the evaluation of bosonic energies, connection probabilities, and related estimators.
 *
 * @see BosonicExchangeBase
 */
class FictitiousExchange final : public BosonicExchangeBase {
public:
    /**
     * @brief Construct a new FictitiousExchange object.
     */
    explicit FictitiousExchange(
        const std::shared_ptr<const VecArray>& coord_first_bead,
        const std::shared_ptr<const VecArray>& coord_last_bead,
        const ThermalContext& thermal_ctx,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        const BeadContext& bead_ctx,
        double exchange_xi
    );

    /**
     * @brief Destructor.
     */
    ~FictitiousExchange() override = default;

    /**
     * @brief Compute the effective spring potential of the entire bosonic system (kinetic part of the density matrix).
     * @return The effective bosonic potential V^[1,N].
     */
    double effectivePotential() override;

    /**
     * @brief Get the effective bosonic potential V^[1,n] for a given number of bosons.
     * @param n Number of bosons ("n" is less than or equal to "N").
     * @return The effective bosonic potential V^[1,n].
     */
    [[nodiscard]] double getVn(int n) const;

    /**
     * @brief Get the cycle energy of the last i particles.
     * @param i Number of particles in the cycle.
     * @return The energy of the given cycle.
     */
    [[nodiscard]] double getEknSerialOrder(int i) const;

    /**
     * @brief Re-evaluate bosonic energies and connection probabilities.
     *        Should be called after coordinate updates.
     */
    void prepare() override;

    /**
     * @brief Evaluates the partial derivative of (beta*V_B) with respect to beta.
     *        This is used to calculate the primitive kinetic energy estimator for bosons.
     *        Corresponds to Eqns. (4)-(5) in SI of pnas.1913365116,
     *        excluding the constant factor of d*N*P/(2*beta).
     * @return The exterior spring contribution to the overall kinetic energy of the quantum system.
     */
    double primitiveEnergyEstimator() override;

    /**
     * @brief Probability of the configuration where all particles are separate (identity permutation).
     * @return Probability of the distinct configuration.
     */
    double getDistinctProbability() override;

    /**
     * @brief Probability of the configuration where all particles are connected (single cycle), divided by 1/N.
     *        Notice that there are (N-1)! permutations of this topology (all represented by the cycle 0,1,...,N-1,0);
     *        this cancels the division by 1/N.
     * @return Probability of the longest cycle configuration.
     */
    double getLongestProbability() override;

    /**
     * @brief The average permutation sign as defined in Eq. (9) https://doi.org/10.1063/5.0008720.
     * @return The average sign.
     */
    double getSign() override;

    /**
     * @brief Print debug information about the bosonic exchange state.
     */
    void printBosonicDebug() override;

protected:
    /**
     * @brief Compute the spring force on the first bead.
     * @param f Output force vector.
     */
    void springForceFirstBead(VecArray& f) override;

    /**
     * @brief Compute the spring force on the last bead.
     * @param f Output force vector.
     */
    void springForceLastBead(VecArray& f) override;

private:
    /**
     * @brief Evaluate all bosonic energies (cycle, prefix, suffix), as well as the connection probabilities.
     */
    void evaluateBosonicEnergies();

    /**
     * @brief Evaluate all the cycle energies, as outlined in Eqs. 5-7 of arXiv:2305.18025.
     */
    void evaluateCycleEnergies();

    /**
     * @brief Get the spring energy of the cycle (m-k+1,...,m)
     * @param m Number of particles.
     * @param k Number of elements in the cycle.
     * @return The energy value.
     */
    [[nodiscard]] double getEnk(int m, int k) const;

    /**
     * @brief Set the spring energy of the cycle (m-k+1,...,m)
     * @param m Number of particles.
     * @param k Number of elements in the cycle.
     * @param val The energy value to set.
     * @return The energy value.
     */
    void setEnk(int m, int k, double val);

    /**
     * @brief Evaluate connection probabilities.
     */
    void evaluateConnectionProbabilities();

    /**
     * @brief Evaluate the prefix (forward) bosonic potentials V^[1,v], as outlined in Eq. 3 of arXiv:2305.18025.
     *        Assumes that the cycle energies have been computed.
     */
    void evaluatePrefixPotential();

    /**
     * @brief Evaluate the suffix (backward) bosonic potentials V^[u,N], as outlined in Eq. 15 of arXiv:2305.18025.
     *        Assumes that both the cycle energies and prefix potentials have been computed.
     */
    void evaluateSuffixPotential();

    double m_xi;                                     ///< The fictitious particle parameter xi
    std::vector<double> m_cycle_energies;            ///< Array of cycle energies E(k,N) = E^[N-k+1,N]
    std::vector<double> m_prefix_pot;                ///< Forward potentials V^[1,1], V^[1,2], ..., V^[1,N]
    std::vector<double> m_suffix_pot;                ///< Backward potentials V^[1,N], V^[2,N], ..., V^[N,N]
    std::vector<double> m_connection_probabilities;  ///< Connection probabilities
    double m_log_n_factorial;                        ///< Logarithm of N! for normalization
};