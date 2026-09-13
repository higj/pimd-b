#pragma once

#include "observables/observable.h"
#include "contexts/bead_context.h"
#include "common.h"
#include <memory>


class CenterOfMassObservable final : public Observable {
public:
    CenterOfMassObservable(
        const std::shared_ptr<const VecArray>& coord,
        const BeadContext& bead_ctx,
        int out_freq,
        const std::string& out_unit
    );

    ~CenterOfMassObservable() override = default;

    void calculate() override;

private:
    std::shared_ptr<const VecArray> m_coord;
    BeadContext m_bead_ctx;
};
