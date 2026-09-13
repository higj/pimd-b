#include "observables/temperature_observable.h"
#include "observables/shared_computations.h"
#include "common.h"

TemperatureObservable::TemperatureObservable(
    const VelocityContext& vel_ctx,
    const BeadContext& bead_ctx,
    int   out_freq,
    const std::string& out_unit
) : Observable("temperature", out_freq, out_unit),
    m_vel_ctx(vel_ctx),
    m_bead_ctx(bead_ctx) 
{
    initializeLabel("temperature");
}

void TemperatureObservable::calculate() {
    // Temperature is calculated according to Tolman's equipartition theorem as the average kinetic
    // energy per degree of freedom. This might not be accurate for systems with constraints.
    // See [J. Chem. Theory Comput. 2019, 15, 1, 84-94.] for a discussion on the topic.
    const double ke = SharedComputations::classicalKineticEnergy(m_vel_ctx, m_bead_ctx, m_cache);

    /// TODO: When zeroing the center of mass motion, the number of degrees of freedom must be reduced by NDIM.
    const double dof = static_cast<double>(NDIM) * m_bead_ctx.natoms * m_bead_ctx.nbeads;

    double temperature = 2.0 * ke / (dof * Constants::kB);

    // In the tau convention, the ring-polymer simulation is performed at P times
    // the physical temperature; divide by nbeads to recover the quantum temperature.
#if TAU_CONVENTION
    temperature /= static_cast<double>(m_bead_ctx.nbeads);
#endif

    quantities["temperature"] = Units::convertToUser("temperature", m_out_unit, temperature);
}