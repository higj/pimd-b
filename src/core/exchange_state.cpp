#include "core/exchange_state.h"
#include "bosonic_exchange.h"

ExchangeState::ExchangeState(
    const std::shared_ptr<const dVec>& coord_first_bead,
    const std::shared_ptr<const dVec>& coord_last_bead,
    const ThermalContext& thermal_ctx,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    int nbeads,
    int this_bead,
    bool bosonic
) :
    is_bosonic(bosonic),
    is_first_bead(this_bead == 0),
    is_last_bead(this_bead == nbeads - 1),
    is_bosonic_bead(bosonic && (this_bead == 0 || this_bead == nbeads - 1))
{
    // If the imaginary time-slice is either 1 or P, initialize the bosonic exchange algorithm
    // CR: otherwise it remains uninitialized?
    // JH: Yes
    if (is_bosonic_bead)
    {
#if FACTORIAL_BOSONIC_ALGORITHM
        bosonic_exchange = std::make_unique<FactorialBosonicExchange>(
            coord_first_bead,
            coord_last_bead,
            thermal_ctx,
            spring_ctx,
            box_ctx,
            nbeads,
            this_bead
        );
#else
        bosonic_exchange = std::make_unique<BosonicExchange>(
            coord_first_bead,
            coord_last_bead,
            thermal_ctx,
            spring_ctx,
            box_ctx,
            nbeads,
            this_bead
        );
#endif
    }
}