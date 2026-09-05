#pragma once

#include "observables/observable.h"
#include "contexts/bead_context.h"
#include "contexts/box_context.h"
#include "contexts/thermal_context.h"
#include "contexts/spring_context.h"
#include "common.h"
#include <memory>

class ForceManager;

class SuzukiChinWeightObservable final : public Observable
{
public:
    SuzukiChinWeightObservable(
        const std::shared_ptr<const VecArray>& coord,
        const std::shared_ptr<const ForceManager>& force_mgr,
        const BeadContext& bead_ctx,
        const BoxContext& box_ctx,
        const ThermalContext& thermal_ctx,
        const SpringContext& spring_ctx,
        double alpha,
        int out_freq,
        const std::string& out_unit
    );

    void calculate() override;

private:
    std::shared_ptr<const VecArray> m_coord;
    std::shared_ptr<const ForceManager> m_force_mgr;
    BeadContext m_bead_ctx;
    BoxContext m_box_ctx;
    ThermalContext m_thermal_ctx;
    SpringContext m_spring_ctx;
    double m_alpha;
};
