#include "core/force_manager.h"
#include "core/statistics_manager.h"
#include "strategies/forces/spring_force_strategy.h"

ForceManager::ForceManager(
    std::unique_ptr<Potential> ext_pot,
    std::unique_ptr<Potential> int_pot,
    double cutoff_in,
    const BoxContext& box_ctx
) : m_spring_force_strategy(StatisticsManager::getInstance().createSpringForceStrategy())
{
    ext_potential = std::move(ext_pot);
    int_potential = std::move(int_pot);
    cutoff = cutoff_in;

    // For cubic cells with PBC, the cutoff must be no greater than L/2
    // to be consistent with the minimum image convention
    // (see section 1.6.3 in Allen & Tildesley).
    if (box_ctx.pbc && cutoff > 0.0)
        cutoff = std::min(cutoff, 0.5 * box_ctx.box_size);
}


ForceManager::~ForceManager() = default;

void ForceManager::updatePhysicalForces(SystemState& state, const BoxContext& box_ctx) const
{
    // Calculate the external forces acting on the particles
    state.physical_forces = (-1.0) * ext_potential->gradV(state.coord);

    if (cutoff == 0.0)
        return;

    const int natoms = state.getNumAtoms();

    for (int ptcl_one = 0; ptcl_one < natoms; ++ptcl_one)
    {
        for (int ptcl_two = ptcl_one + 1; ptcl_two < natoms; ++ptcl_two)
        {
            // Get the vector distance between the two particles.
            Vec diff = state.coord.getSeparationArray(ptcl_one, ptcl_two);
            /// TODO: MINIM should become a parameter (mic_spring and mic_potential)
            box_ctx.applyMinimumImageIfNeeded(diff);

            // If the distance between the particles exceeds the cutoff length
            // then we assume the interaction is negligible and do not bother
            // calculating the force.
            // We use the convention that when cutoff < 0 then the interaction is
            // calculated for all distances.
            if (const double distance = norm(diff); distance < cutoff || cutoff < 0.0)
            {
                Vec grad_on_one = int_potential->gradV(diff);

                for (int axis = 0; axis < NDIM; ++axis)
                {
                    state.physical_forces(ptcl_one, axis) -= grad_on_one[axis];
                    state.physical_forces(ptcl_two, axis) += grad_on_one[axis];
                }
            }
        }
    }
}

void ForceManager::updateForces(SystemState& state, const SpringContext& spring_ctx, const BoxContext& box_ctx) const
{
    // First, update the spring forces based on the current state of the system.
    m_spring_force_strategy->updateSpringForces(state, spring_ctx, box_ctx);

    // Then, update the physical forces acting on the particles.
    updatePhysicalForces(state, box_ctx);
}
