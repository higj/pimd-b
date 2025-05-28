#include <array>

#include "bosonic_exchange/bosonic_exchange_base.h"

BosonicExchangeBase::BosonicExchangeBase(const BosonicExchangeContext& context) : m_context(context)
{
}

void BosonicExchangeBase::getExteriorBeadsSeparation(int first_idx, int last_idx, std::array<double, NDIM>& diff) const
{
    first_idx = first_idx % m_context.nbosons;
    last_idx = last_idx % m_context.nbosons;

    for (int axis = 0; axis < NDIM; ++axis)
    {
        double dx = (*m_context.x_first_bead)(first_idx, axis) - (*m_context.x_last_bead)(last_idx, axis);
#if MINIM
        if (m_context.pbc)
            applyMinimumImage(dx, m_context.box_size);
#endif
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
void BosonicExchangeBase::exteriorSpringForce(dVec& f)
{
    if (m_context.this_bead == 0)
    {
        springForceFirstBead(f);
    }
    else
    {
        springForceLastBead(f);
    }
}
