#pragma once

#include "initializable_item.h"

#include <cstdint>

enum class DumpType : std::uint8_t {
    POSITION,
    VELOCITY,
    FORCE,
    UNKNOWN
};

class DumpItem final : public InitializableItem<DumpType> {
public:
    DumpItem(const std::string& dump_name, const std::string& dump_unit) {
        name = dump_name;
        unit = dump_unit;
    }

    [[nodiscard]] DumpType getType() const override {
        if (name == "positions") return DumpType::POSITION;
        if (name == "velocities") return DumpType::VELOCITY;
        if (name == "forces") return DumpType::FORCE;
        return DumpType::UNKNOWN;
    }
};
