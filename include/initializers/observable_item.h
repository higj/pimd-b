#pragma once

#include "initializable_item.h"

#include <cstdint>

enum class ObservableType : std::uint8_t {
    ENERGY,
    CLASSICAL,
    BOSONIC,
    GSF,
    UNKNOWN
};

class ObservableItem final : public InitializableItem<ObservableType> {
public:
    ObservableItem(const std::string& obs_name, const std::string& obs_unit) {
        name = obs_name;
        unit = obs_unit;
    }

    [[nodiscard]] ObservableType getType() const override {
        if (name == "energy") return ObservableType::ENERGY;
        if (name == "classical") return ObservableType::CLASSICAL;
        if (name == "bosonic") return ObservableType::BOSONIC;
        if (name == "gsf") return ObservableType::GSF;
        return ObservableType::UNKNOWN;
    }
};
