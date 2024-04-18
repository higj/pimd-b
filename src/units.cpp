#include <stdexcept>
#include <regex>
#include <ranges>

#include "units.h"


namespace Units {
    std::tuple<std::string, std::string> separatePrefixUnit(const std::string& unit) {
        std::string regex_pattern;

        for (const std::string& key : UnitPrefix | std::views::keys) {
            if (!key.empty()) {
                regex_pattern += key + "|";
            }
        }

        regex_pattern.pop_back(); // Remove last "|"
        regex_pattern = "(" + regex_pattern + ")*(.*)";

        static const std::regex prefix_regex(regex_pattern);

        if (std::smatch match; std::regex_match(unit, match, prefix_regex)) {
            // match[0] represents the full string
            // match[1] is the prefix
            // match[2] is the unit
            return std::make_tuple(match[1].str(), match[2].str());
        }

        // No prefix found, return empty string for prefix and the input unit
        return std::make_tuple("", unit);
    }

    double convertToInternal(const std::string& family, const std::string& unit, const double number) {
        if (family == "number") {
            return number;
        }

        auto it_family = UnitMap.find(family);
        if (it_family == UnitMap.end()) {
            throw std::invalid_argument(family + " is an undefined units kind.");
        }

        auto [prefix, base_unit] = separatePrefixUnit(unit);

        auto it_base = it_family->second.find(base_unit);
        if (it_base == it_family->second.end()) {
            throw std::invalid_argument(base_unit + " is an undefined unit for kind " + family + ".");
        }

        return number * it_base->second * UnitPrefix.at(prefix);
    }

    double convertToUser(const std::string& family, const std::string& unit, double number) {
        return number / convertToInternal(family, unit, 1.0);
    }
}
