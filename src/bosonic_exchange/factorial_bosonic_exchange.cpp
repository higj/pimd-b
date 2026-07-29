#include <array>
#include <numeric>
#include <algorithm>

#include "bosonic_exchange/factorial_bosonic_exchange.h"

FactorialBosonicExchange::FactorialBosonicExchange(
    const std::shared_ptr<const VecArray>& coord_first_bead,
    const std::shared_ptr<const VecArray>& coord_last_bead,
    const ThermalContext& thermal_ctx,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    const BeadContext& bead_ctx
) : BosonicExchangeBase(
        coord_first_bead,
        coord_last_bead,
        thermal_ctx,
        spring_ctx,
        box_ctx,
        bead_ctx
    ),
    m_labels(bead_ctx.natoms)
{
    // Fill the labels array with numbers from 0 to nbosons-1
    std::iota(m_labels.begin(), m_labels.end(), 0);

    // For numerical stability
    m_e_shift = getMinExteriorSpringEnergy();
}

void FactorialBosonicExchange::prepare()
{
    m_e_shift = getMinExteriorSpringEnergy();
}

int FactorialBosonicExchange::firstBeadNeighbor(int ptcl_idx) const
{
    return static_cast<int>(std::distance(m_labels.begin(), std::ranges::find(m_labels, ptcl_idx)));
}

int FactorialBosonicExchange::lastBeadNeighbor(int ptcl_idx) const
{
    return m_labels[ptcl_idx];
}

double FactorialBosonicExchange::getMinExteriorSpringEnergy()
{
    double min_delta = std::numeric_limits<double>::max();

    // Iterate over all permutations
    do
    {
        double diff2 = 0.0;

        for (int l = 0; l < m_bead_ctx.natoms; ++l)
        {
            // First bead of some particle (depending on the permutation) minus last bead of the l-th particle
            diff2 += getExteriorSeparationSquared(l, lastBeadNeighbor(l));
        }

        // Compare the current total exterior spring energy with the minimum total exterior spring energy found so far
        min_delta = std::min(min_delta, 0.5 * m_spring_ctx.spring_constant * diff2);
    }
    while (std::ranges::next_permutation(m_labels).found);

    return min_delta;
}

double FactorialBosonicExchange::effectivePotential()
{
    long permutation_counter = 0;
    double weights_sum = 0.0;

    // Iterate over all permutations and calculate the weights
    // associated with the exterior beads.
    do
    {
        permutation_counter++;
        double diff2 = 0.0;

        for (int ptcl_idx = 0; ptcl_idx < m_bead_ctx.natoms; ++ptcl_idx)
        {
            diff2 += getExteriorSeparationSquared(ptcl_idx, lastBeadNeighbor(ptcl_idx));
        }

        weights_sum += exp(-m_spring_ctx.beta_half_k * diff2);
    }
    while (std::ranges::next_permutation(m_labels).found);

    return (-1.0 / m_thermal_ctx.thermo_beta) * log(weights_sum / permutation_counter);
}

void FactorialBosonicExchange::springForceLastBead(VecArray& f)
{
    /// TODO: Either reset "f" at the beginning of each MD step, or don't "+=" the force later in this function
    f.reset();

    VecArray temp_force(m_bead_ctx.natoms);
    double denom_weight = 0.0;

    do
    {
        double weight = 0.0;

        for (int l = 0; l < m_bead_ctx.natoms; ++l)
        {
            std::array<double, NDIM> diff_next;

            // First bead (1) of some particle (depending on the permutation) minus last bead (P) of the l-th particle
            getExteriorBeadsSeparation(lastBeadNeighbor(l), l, diff_next);

            // Coordinate differences corresponding to exterior beads
            for (int axis = 0; axis < NDIM; ++axis)
            {
                weight += diff_next[axis] * diff_next[axis];
                temp_force(l, axis) = diff_next[axis] * m_spring_ctx.spring_constant;
            }
        }

        weight = exp(-m_thermal_ctx.thermo_beta * (0.5 * m_spring_ctx.spring_constant * weight - m_e_shift));

        for (int l = 0; l < m_bead_ctx.natoms; ++l)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                f(l, axis) += weight * temp_force(l, axis);
            }
        }

        denom_weight += weight;
    }
    while (std::ranges::next_permutation(m_labels).found);

    // Normalize the forces by the sum of Boltzmann contributions due to all the permutations
    for (int l = 0; l < m_bead_ctx.natoms; ++l)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            f(l, axis) = f(l, axis) / denom_weight;
        }
    }
}

void FactorialBosonicExchange::springForceFirstBead(VecArray& f)
{
    /// TODO: Either reset "f" at the beginning of each MD step, or don't "+=" the force later in this function
    f.reset();

    VecArray temp_force(m_bead_ctx.natoms);
    double denom_weight = 0.0;

    do
    {
        double weight = 0.0;

        for (int l = 0; l < m_bead_ctx.natoms; ++l)
        {
            std::array<double, NDIM> diff_prev;

            // Last bead (P) of some particle (depending on the permutation) minus first bead (1) of the l-th particle
            getExteriorBeadsSeparation(l, firstBeadNeighbor(l), diff_prev);

            // Coordinate differences corresponding to exterior beads
            for (int axis = 0; axis < NDIM; ++axis)
            {
                weight += diff_prev[axis] * diff_prev[axis];
                temp_force(l, axis) = -diff_prev[axis] * m_spring_ctx.spring_constant;
            }
        }

        weight = exp(-m_thermal_ctx.thermo_beta * (0.5 * m_spring_ctx.spring_constant * weight - m_e_shift));

        for (int l = 0; l < m_bead_ctx.natoms; ++l)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                f(l, axis) += weight * temp_force(l, axis);
            }
        }

        denom_weight += weight;
    }
    while (std::ranges::next_permutation(m_labels).found);

    // Normalize the forces by the sum of Boltzmann contributions due to all the permutations
    for (int l = 0; l < m_bead_ctx.natoms; ++l)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            f(l, axis) = f(l, axis) / denom_weight;
        }
    }
}

double FactorialBosonicExchange::getDistinctProbability()
{
    /// @todo Currently not implemented
    return 0.0;
}

double FactorialBosonicExchange::getLongestProbability()
{
    /// @todo Currently not implemented
    return 0.0;
}

double FactorialBosonicExchange::getSign()
{
    /// @todo Currently not implemented
    return 0.0;
}

double FactorialBosonicExchange::primitiveEnergyEstimator()
{
    double numerator = 0.0;
    double denom_weight = 0.0;

    do
    {
        double weight = 0.0;

        for (int l = 0; l < m_bead_ctx.natoms; ++l)
        {
            // Last bead of some particle (depending on the permutation) minus first bead of the l-th particle
            weight += getExteriorSeparationSquared(l, firstBeadNeighbor(l));
        }

        const double delta_spring_energy = 0.5 * m_spring_ctx.spring_constant * weight;
        weight = exp(-m_thermal_ctx.thermo_beta * (delta_spring_energy - m_e_shift));

        // Delta(E^sigma) is ready only at this point. We multiply Delta(E^sigma) by the weight
        // and add the result to the numerator (and also update the denominator by adding the 
        // weight of the current permutation).
        numerator += delta_spring_energy * weight;
        denom_weight += weight;
    }
    while (std::ranges::next_permutation(m_labels).found);

    double bosonic_spring_energy = numerator / denom_weight;

#if IPI_CONVENTION
    bosonic_spring_energy /= m_bead_ctx.nbeads;
#endif

    return (-1.0) * bosonic_spring_energy;
}

void FactorialBosonicExchange::printBosonicDebug()
{
}
