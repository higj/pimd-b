#pragma once
#include "observables/observable.h"
#include "contexts/bead_context.h"
#include "common.h"

#include <memory>

class VirialEnergyObservable final : public Observable {
public:
    VirialEnergyObservable(
        const std::shared_ptr<const VecArray>& coord,
        const std::shared_ptr<const VecArray>& physical_forces,
        const BeadContext& bead_ctx,
        int out_freq,
        const std::string& out_unit
    );

    void calculate() override;

private:
    std::shared_ptr<const VecArray> m_coord;
    std::shared_ptr<const VecArray> m_forces;
    BeadContext m_bead_ctx;
};
