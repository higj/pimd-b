#include "strategies/observables/distinguishable_classical_spring_energy_strategy.h"
#include "ring_polymer_utils.h"

double DistinguishableClassicalSpringEnergyStrategy::calculateSpringEnergy(
    const dVec& coord,
    const dVec& prev_coord,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx
)
{
    return RingPolymerUtils::classicalSpringEnergy(
        coord,
        prev_coord,
        spring_ctx.spring_constant,
        box_ctx
    );
}
