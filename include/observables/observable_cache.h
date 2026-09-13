#pragma once

#include "observables/cache_keys.h"

#include <array>
#include <optional>

/**
 * @brief Step-scoped cache for intermediate observable quantities.
 *
 * Owned by ObservablesLogger. Invalidated once per logging step, before
 * any observable's calculate() is called. Observables that share an
 * intermediate computation write to and read from this cache using 
 * agreed-upon enum values.
 *
 * All quantities are stored in internal (unconverted) units. Unit conversion
 * happens in the observable after the value is retrieved.
 */
class ObservableCache {
public:
    void invalidate() { m_values.fill(std::nullopt); }

    [[nodiscard]] std::optional<double> get(CacheKey key) const {
        return m_values[static_cast<std::size_t>(key)];
    }

    void put(CacheKey key, double value) {
        m_values[static_cast<std::size_t>(key)] = value;
    }

private:
    std::array<std::optional<double>,
        static_cast<std::size_t>(CacheKey::COUNT)> m_values{};
};