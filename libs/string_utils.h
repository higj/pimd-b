#pragma once

#include <string>
#include <vector>

// Define a list of strings
using StringsList = std::vector<std::string>;

namespace StringUtils {
    bool parseTokenParentheses(const std::string& input, std::string& token, std::string& value);
    StringsList splitString(const std::string& line, char delimiter = ' ');
    bool labelInArray(const std::string& label, const StringsList& arr);
}
