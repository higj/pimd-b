#pragma once

#include "bosonic_exchange/bosonic_exchange_base.h"

#include <vector>

class FactorialBosonicExchange final : public BosonicExchangeBase {
public:
    explicit FactorialBosonicExchange(
        const std::shared_ptr<const VecArray>& coord_first_bead,
        const std::shared_ptr<const VecArray>& coord_last_bead,
        const ThermalContext& thermal_ctx,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        const BeadContext& bead_ctx
    );
    ~FactorialBosonicExchange() override = default;

    /**
     * @brief Recalculates the smallest exterior spring energy after each coordinate update.
     */
    void prepare() override;

    /**
     * Calculates the effective bosonic spring potential. This is V_eff from the exp(-beta * V_eff) in the
     * partition function. In contrast to the classical spring energy of interior beads,
     * this contribution is no longer harmonic, but is a sum of Boltzmann weights due to all permutations.
     *
     * @return Effective bosonic exchange potential.
     */
    double effectivePotential() override;

    /**
     * Evaluate the probability of the configuration where all the particles are separate.
     *
     * @return Probability of a configuration corresponding to the identity permutation.
     */
    double getDistinctProbability() override;

    /**
     * Evaluate the probability of a configuration where all the particles are connected,
     * divided by 1/N. Notice that there are (N-1)! permutations of this topology
     * (all represented by the cycle 0,1,...,N-1,0); this cancels the division by 1/N.
     *
     * @return Probability of a configuration where all the particles are connected.
     */
    double getLongestProbability() override;

    /**
     * Calculates the contribution of the exterior beads to the primitive kinetic energy estimator.
     *
     * @return Weighted average of exterior spring energies over all permutations.
     */
    double primitiveEnergyEstimator() override;

    void printBosonicDebug() override;

protected:
    /**
     * Evaluates the bosonic spring forces acting on the particles at the last imaginary time slice.
     *
     * @param[out] f Spring forces acting on the particles at time-slice P.
     */
    void springForceLastBead(VecArray& f) override;

    /**
     * Evaluates the bosonic spring forces acting on the particles at the first imaginary time slice.
     *
     * @param[out] f Spring forces acting on the particles at time-slice 1.
     */
    void springForceFirstBead(VecArray& f) override;

private:
    /**
     * Identifies the particle to which the neighboring P bead, connected to the first bead of ptcl_idx, belongs.
     * In the context of Cauchy's two-line notation, it retrieves the number from the upper line, placed above ptcl_idx.
     * Equivalently, it is the position of ptcl_idx in the bottom line, assuming zero-based indexing.
     * This works because the upper line can be viewed as the list of particles in the last imaginary time slice and
     * the bottom line as the list of their respective neighbors in the first imaginary time slice.
     *
     * @param ptcl_idx Index of the particle associated with the first bead.
     * @return Index of the particle associated with the neighboring P bead.
     */
    [[nodiscard]] int firstBeadNeighbor(int ptcl_idx) const;

    /**
     * Identifies the particle to which the neighboring 1 bead, connected to the last bead of ptcl_idx, belongs.
     * In the context of Cauchy's two-line notation, it retrieves the number from the bottom line, placed below ptcl_idx.
     * Equivalently, it is the particle label in the bottom line at position ptcl_idx, assuming zero-based indexing.
     * This works because the upper line can be viewed as the list of particles in the last imaginary time slice and
     * the bottom line as the list of their respective neighbors in the first imaginary time slice.
     *
     * @param ptcl_idx Index of the particle associated with the last bead.
     * @return Index of the particle associated with the neighboring 1 bead.
     */
    [[nodiscard]] int lastBeadNeighbor(int ptcl_idx) const;

    /**
     * Calculates the smallest exterior spring energy that a permutation can yield in a given time-step. This energy is then used to shift
     * the spring energies in the Boltzmann weights to avoid numerical instabilities, which can be especially prominent
     * at high Trotter numbers.
     *
     * @return The smallest exterior spring energy contribution due to a permutation.
     */
    double getMinExteriorSpringEnergy();

    std::vector<int> m_labels;  ///< Particle labels
    double m_e_shift;           ///< Shift of the exterior spring energy, used to avoid numerical instabilities
};
