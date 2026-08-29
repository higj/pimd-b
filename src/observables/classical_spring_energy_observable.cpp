#include "observables/classical_spring_energy_observable.h"
#include "core/statistics_manager.h"
#include "strategies/observables/classical_spring_energy_strategy.h"

ClassicalSpringEnergyObservable::ClassicalSpringEnergyObservable(
    const std::shared_ptr<const VecArray>& coord,
    const std::shared_ptr<const VecArray>& prev_coord,
    const SpringContext& spring_ctx,
    const BoxContext& box_ctx,
    int   out_freq,
    const std::string& out_unit
) : Observable("cl_spring", out_freq, out_unit),
    m_spring_energy_strategy(StatisticsManager::getInstance().createClassicalSpringEnergyStrategy()),
    m_coord(coord),
    m_prev_coord(prev_coord),
    m_spring_ctx(spring_ctx),
    m_box_ctx(box_ctx)
{
    initializeLabel("cl_spring");
}

ClassicalSpringEnergyObservable::~ClassicalSpringEnergyObservable() = default;

void ClassicalSpringEnergyObservable::calculate() {
    const double spring_energy = m_spring_energy_strategy->calculateSpringEnergy(
        *m_coord,
        *m_prev_coord,
        m_spring_ctx,
        m_box_ctx
    );
    quantities["cl_spring"] = Units::convertToUser("energy", m_out_unit, spring_energy);
}