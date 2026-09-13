#include "observables/primitive_kinetic_energy_observable.h"
#include "core/statistics_manager.h"
#include "strategies/observables/primitive_kinetic_energy_strategy.h"
#include "common.h"

PrimitiveKineticEnergyObservable::PrimitiveKineticEnergyObservable(
    const std::shared_ptr<const VecArray>& coord,
    const std::shared_ptr<const VecArray>& prev_coord,
    const BeadContext& bead_ctx,
    const ThermalContext& thermal_ctx,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    int out_freq,
    const std::string& out_unit
) : Observable("kinetic", out_freq, out_unit),
m_strategy(StatisticsManager::getInstance().createPrimitiveKineticEnergyStrategy()),
m_coord(coord),
m_prev_coord(prev_coord),
m_bead_ctx(bead_ctx),
m_thermal_ctx(thermal_ctx),
m_spring_ctx(spring_ctx),
m_box_ctx(box_ctx) {
    initializeLabel("kinetic");
}

PrimitiveKineticEnergyObservable::~PrimitiveKineticEnergyObservable() = default;

void PrimitiveKineticEnergyObservable::calculate() {
    double ke = 0.5 * NDIM * m_bead_ctx.natoms / m_thermal_ctx.beta;
    ke += m_strategy->calculateSpringContribution(
        *m_coord, *m_prev_coord, m_spring_ctx, m_box_ctx, m_bead_ctx
    );
    quantities["kinetic"] = Units::convertToUser("energy", m_out_unit, ke);
}