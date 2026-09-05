#pragma once
#include "observables/observable.h"
#include "contexts/bead_context.h"
#include "contexts/thermal_context.h"
#include "contexts/spring_context.h"
#include "contexts/box_context.h"
#include <memory>

class PrimitiveKineticEnergyStrategy;

class PrimitiveKineticEnergyObservable final : public Observable {
public:
    PrimitiveKineticEnergyObservable(
        const std::shared_ptr<const VecArray>& coord,
        const std::shared_ptr<const VecArray>& prev_coord,
        const BeadContext& bead_ctx,
        const ThermalContext& thermal_ctx,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        int out_freq,
        const std::string& out_unit
    );
    ~PrimitiveKineticEnergyObservable() override;

    void calculate() override;

private:
    std::unique_ptr<PrimitiveKineticEnergyStrategy> m_strategy;
    std::shared_ptr<const VecArray> m_coord;
    std::shared_ptr<const VecArray> m_prev_coord;
    BeadContext    m_bead_ctx;
    ThermalContext m_thermal_ctx;
    SpringContext  m_spring_ctx;
    BoxContext     m_box_ctx;
};