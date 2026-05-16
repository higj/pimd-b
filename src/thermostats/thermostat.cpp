#include "thermostats/thermostat.h"
#include "core/system_state.h"
#include "thermostats/thermostat_coupling.h"

Thermostat::Thermostat(
    const ThermalContext& thermal_ctx, 
    const NormalModesContext& nm_ctx,
    const std::shared_ptr<SystemState>& state
) : m_thermal_ctx(thermal_ctx), m_nm_ctx(nm_ctx), m_state(state)
{
    const auto& momenta_ptr = std::shared_ptr<VecArray>(state, &state->momenta);

    // Choose coupling (Cartesian coords or normal modes of distinguishable ring polymers)
    if (nm_ctx.couple_to_nm)
    {
        m_coupling = std::make_unique<NormalModesCoupling>(
            momenta_ptr,
            nm_ctx.normal_modes,
            state->currentBead()
        );
    }
    else
    {
        m_coupling = std::make_unique<CartesianCoupling>(
            momenta_ptr
        );
    }
}

Thermostat::~Thermostat() = default;

// This is the step function of a general thermostat, called in the simulation's run loop
void Thermostat::step()
{
    m_coupling->mpiCommunication();
    momentaUpdate();
    m_coupling->updateCoupledMomenta();
}

// This is an update of the momenta within the thermostat step, unique for each thermostat
void Thermostat::momentaUpdate()
{
}

double Thermostat::getAdditionToH()
{
    return 0;
}
