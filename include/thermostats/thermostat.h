#pragma once

#include "contexts/thermal_context.h"
#include "contexts/normal_modes_context.h"

class SystemState;
class Coupling;

class Thermostat
{
public:
    explicit Thermostat(
        const ThermalContext& thermal_ctx,
        const NormalModesContext& nm_ctx,
        const std::shared_ptr<SystemState>& state
    );
    virtual ~Thermostat() = default;

    void step();
    virtual void momentaUpdate();
    virtual double getAdditionToH();

protected:
    ThermalContext m_thermal_ctx;
    NormalModesContext m_nm_ctx;
    std::shared_ptr<SystemState> m_state;
    std::unique_ptr<Coupling> m_coupling;
};
