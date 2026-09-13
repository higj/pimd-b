#include "observables/classical_kinetic_energy_observable.h"
#include "observables/shared_computations.h"

ClassicalKineticEnergyObservable::ClassicalKineticEnergyObservable(
    const VelocityContext& vel_ctx,
    const BeadContext& bead_ctx,
    int   out_freq,
    const std::string& out_unit
) : Observable("cl_kinetic", out_freq, out_unit),
    m_vel_ctx(vel_ctx),
    m_bead_ctx(bead_ctx) 
{
    initializeLabel("cl_kinetic");
}

void ClassicalKineticEnergyObservable::calculate() {
    const double ke = SharedComputations::classicalKineticEnergy(m_vel_ctx, m_bead_ctx, m_cache);
    quantities["cl_kinetic"] = Units::convertToUser("energy", m_out_unit, ke);
}