#include <array>

#include "bosonic_exchange/bosonic_exchange_base.h"

BosonicExchangeBase::BosonicExchangeBase(
    const std::shared_ptr<const VecArray>& coord_first_bead,
    const std::shared_ptr<const VecArray>& coord_last_bead,
    const ThermalContext& thermal_ctx,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    const BeadContext& bead_ctx,
    double exchange_xi
) : m_coord_first_bead(coord_first_bead),
    m_coord_last_bead(coord_last_bead),
    m_thermal_ctx(thermal_ctx),
    m_spring_ctx(spring_ctx),
    m_box_ctx(box_ctx),
    m_bead_ctx(bead_ctx),
    m_exchange_xi(exchange_xi)
{
}

void BosonicExchangeBase::getExteriorBeadsSeparation(int first_idx, int last_idx, std::array<double, NDIM>& diff) const
{
    first_idx = first_idx % m_bead_ctx.natoms;
    last_idx = last_idx % m_bead_ctx.natoms;

    for (int axis = 0; axis < NDIM; ++axis)
    {
        double dx = (*m_coord_first_bead)(first_idx, axis) - (*m_coord_last_bead)(last_idx, axis);
        m_box_ctx.applyMinimumImageIfNeeded(dx);

        diff[axis] = dx;
    }
}

double BosonicExchangeBase::getExteriorSeparationSquared(int first_idx, int last_idx) const
{
    std::array<double, NDIM> diff;
    getExteriorBeadsSeparation(first_idx, last_idx, diff);

    double dist_sqrd = 0.0;

    for (int axis = 0; axis < NDIM; ++axis)
    {
        dist_sqrd += diff[axis] * diff[axis];
    }

    return dist_sqrd;
}

/**
 * Calculates the bosonic force exerted on the beads
 * at the first and the last time-slices.
 *
 * @param[out] f Vector to store the forces.
 */
void BosonicExchangeBase::exteriorSpringForce(VecArray& f)
{
    if (m_bead_ctx.this_bead == 0)
    {
        springForceFirstBead(f);
    }
    else
    {
        springForceLastBead(f);
    }
}
