#pragma once

#include "initializable_item.h"

#include <cstdint>

enum class ObservableType : std::uint8_t {
    ENERGY,
    CLASSICAL,
    CL_KINETIC_ENERGY,
    QUANTUM_TEMPERATURE,
    BOSONIC,
    GSF,
    CENTER_OF_MASS,
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
        if (name == "cl_kinetic") return ObservableType::CL_KINETIC_ENERGY;
        if (name == "temperature") return ObservableType::QUANTUM_TEMPERATURE;
        if (name == "bosonic") return ObservableType::BOSONIC;
        if (name == "gsf") return ObservableType::GSF;
        if (name == "center_of_mass") return ObservableType::CENTER_OF_MASS;
        return ObservableType::UNKNOWN;
    }
};
