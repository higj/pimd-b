#include "units.h"

namespace Units {
    std::tuple<std::string, std::string> separate_prefix_unit(const std::string& unit) {
        std::string regex_pattern = "";

        for (auto& it : UnitPrefix) {
            if (it.first != "") {
                regex_pattern += it.first + "|";
            }
        }

        regex_pattern.pop_back(); // Remove last "|"
        regex_pattern = "(" + regex_pattern + ")*(.*)";

        static const std::regex prefixRegex(regex_pattern);

        std::smatch match;
        if (std::regex_match(unit, match, prefixRegex)) {
            // match[0] represents the full string
            // match[1] is the prefix
            // match[2] is the unit
            return std::make_tuple(match[1].str(), match[2].str());
        }
        else {
            // No prefix found, return empty string for prefix and the input unit
            return std::make_tuple("", unit);
        }
    }

    double unit_to_internal(const std::string& family, const std::string& unit, double number) {
        if (family == "number") {
            return number;
        }

        auto itFamily = UnitMap.find(family);
        if (itFamily == UnitMap.end()) {
            throw std::invalid_argument(family + " is an undefined units kind.");
        }

        auto [prefix, base_unit] = separate_prefix_unit(unit);

        auto itBase = itFamily->second.find(base_unit);
        if (itBase == itFamily->second.end()) {
            throw std::invalid_argument(base_unit + " is an undefined unit for kind " + family + ".");
        }

        return number * itBase->second * UnitPrefix.at(prefix);
    }

    double unit_to_user(const std::string& family, const std::string& unit, double number) {
        return number / unit_to_internal(family, unit, 1.0);
    }
}