#pragma once

#include <string>
#include <vector>
#include <array>
#include <string_view>
#include <algorithm>

// Define a list of strings
using StringsList = std::vector<std::string>;

namespace StringUtils {
    bool parseTokenParentheses(const std::string& input, std::string& token, std::string& value);
    StringsList splitString(const std::string& line, char delimiter = ' ');
    bool labelInArray(const std::string& label, const StringsList& arr);

    template <size_t N>
    bool labelInArray(std::string_view label, const std::array<std::string_view, N>& arr) {
        return std::ranges::find(arr, label) != arr.end();
    }
}
