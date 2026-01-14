#pragma once

#include <string>

// Base interface for items that can be enabled/disabled and have physical dimensions.
template<typename ItemType>
class InitializableItem {
public:
    std::string name;
    std::string unit;

    virtual ~InitializableItem() = default;

    [[nodiscard]] bool isEnabled() const
    {
        return (unit != "off" && unit != "false");
    }

    [[nodiscard]] std::string getEffectiveUnit() const
    {
        if (unit == "true" || unit == "on") return "atomic_unit";  // Default unit
        if (unit == "none") return "";                             // No unit conversion
        return unit;                                               // Use specified unit
    }

    [[nodiscard]] virtual ItemType getType() const = 0;
};

/*
#include "initializers/initializable_item.h"

template<typename ItemType>
bool InitializableItem<ItemType>::isEnabled() const {
    return (unit != "off" && unit != "false");
}

template<typename ItemType>
std::string InitializableItem<ItemType>::getEffectiveUnit() const {
    if (unit == "true" || unit == "on") return "atomic_unit";
    if (unit == "none") return "";
    return unit;
}
*/