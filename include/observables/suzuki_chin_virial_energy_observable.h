#pragma once
#include "observables/observable.h"
#include "contexts/bead_context.h"
#include "common.h"
#include <memory>

class ForceManager;

class SuzukiChinVirialEnergyObservable final : public Observable {
public:
    SuzukiChinVirialEnergyObservable(
        const std::shared_ptr<const VecArray>& coord,
        const std::shared_ptr<const ForceManager>& force_mgr,
        const BeadContext& bead_ctx,
        int out_freq,
        const std::string& out_unit
    );

    void calculate() override;

private:
    std::shared_ptr<const VecArray>     m_coord;
    std::shared_ptr<const ForceManager> m_force_mgr;
    BeadContext m_bead_ctx;
};