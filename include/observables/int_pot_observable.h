#pragma once

#include "observables/observable.h"
#include "contexts/bead_context.h"
#include "contexts/box_context.h"
#include "common.h"

#include <memory>

class ForceManager;

class IntPotObservable final : public Observable {
public:
    IntPotObservable(
        const std::shared_ptr<const VecArray>& coord,
        const std::shared_ptr<const ForceManager>& force_mgr,
        const BoxContext& box_ctx,
        const BeadContext& bead_ctx,
        int out_freq,
        const std::string& out_unit
    );

    void calculate() override;

private:
    std::shared_ptr<const VecArray>     m_coord;
    std::shared_ptr<const ForceManager> m_force_mgr;
    BoxContext  m_box_ctx;
    BeadContext m_bead_ctx;
};