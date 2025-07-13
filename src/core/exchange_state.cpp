#include "core/exchange_state.h"
#include "bosonic_exchange.h"

ExchangeState::ExchangeState(
    const std::shared_ptr<const dVec>& coord_first_bead,
    const std::shared_ptr<const dVec>& coord_last_bead,
    const ThermalContext& thermal_ctx,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    const BeadContext& bead_ctx,
    bool bosonic
) :
    is_bosonic(bosonic),
    is_first_bead(bead_ctx.this_bead == 0),
    is_last_bead(bead_ctx.this_bead == bead_ctx.nbeads - 1),
    is_bosonic_bead(bosonic && (bead_ctx.this_bead == 0 || bead_ctx.this_bead == bead_ctx.nbeads - 1))
{
    // If the imaginary time-slice is either 1 or P, initialize the bosonic exchange algorithm
    // CR: otherwise it remains uninitialized?
    // JH: Yes
    // CR: Don't leave fields uninitialized.
    // CR: In any case, this class is a very thin wrapper for BosonicExchange. It can be removed
    if (is_bosonic_bead)
    {
#if FACTORIAL_BOSONIC_ALGORITHM
        bosonic_exchange = std::make_unique<FactorialBosonicExchange>(
            coord_first_bead,
            coord_last_bead,
            thermal_ctx,
            spring_ctx,
            box_ctx,
            bead_ctx
        );
#else
        bosonic_exchange = std::make_unique<BosonicExchange>(
            coord_first_bead,
            coord_last_bead,
            thermal_ctx,
            spring_ctx,
            box_ctx,
            bead_ctx
        );
#endif
    }
}