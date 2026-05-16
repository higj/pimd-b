#include "strategies/observables/distinguishable_primitive_kinetic_energy_strategy.h"
#include "ring_polymer_utils.h"

double DistinguishablePrimitiveKineticEnergyStrategy::calculateSpringContribution(
    const VecArray& coord,
    const VecArray& prev_coord,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    const BeadContext& bead_ctx
)
{
    double spring_energy = RingPolymerUtils::classicalSpringEnergy(
        coord,
        prev_coord,
        spring_ctx.spring_constant,
        box_ctx
    );

#if IPI_CONVENTION
    spring_energy /= bead_ctx.nbeads;
#endif

    return -spring_energy;
}
