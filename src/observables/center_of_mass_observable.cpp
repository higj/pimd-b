#include "observables/center_of_mass_observable.h"

CenterOfMassObservable::CenterOfMassObservable(
    const std::shared_ptr<const VecArray>& coord,
    const BeadContext& bead_ctx,
    int out_freq,
    const std::string& out_unit
) : Observable("center_of_mass", out_freq, out_unit),
    m_coord(coord),
    m_bead_ctx(bead_ctx)
{
    std::vector<std::string> labels;
    labels.reserve(NDIM);

    // Each axis is printed separately
    for (int axis = 0; axis < NDIM; ++axis) {
        labels.push_back("com_" + std::to_string(axis));
    }

    initialize(labels);
}

void CenterOfMassObservable::calculate()
{
    Vec local_com{};
    const auto& coord = *m_coord;

    for (int ptcl_idx = 0; ptcl_idx < m_bead_ctx.natoms; ++ptcl_idx)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            local_com[axis] += coord(ptcl_idx, axis);
        }
    }

    for (int axis = 0; axis < NDIM; ++axis)
    {
        // Normalize by the number of quantum particles
        local_com[axis] /= m_bead_ctx.natoms;

        // Store the COM of this replica
        quantities["com_" + std::to_string(axis)] = Units::convertToUser("length", m_out_unit, local_com[axis]);
    }
}
