#include "string_utils.h"
#include <algorithm>
#include <regex>
#include <sstream>

namespace StringUtils {
    bool labelInArray(const std::string& label, const StringsList& arr) {
        return std::ranges::find(arr, label) != arr.end();
    }

    // Splits a string into several strings based on a delimiter
    std::vector<std::string> splitString(const std::string& line, char delimiter) {
        std::stringstream ss(line);
        std::vector<std::string> tokens;
        std::string token;

        while (std::getline(ss, token, delimiter)) {
            if (!token.empty())
                tokens.push_back(token);
        }

        return tokens;
    }

    // Used to parse strings of the format "token(value)"
    bool parseTokenParentheses(const std::string& input, std::string& token, std::string& value) {
        // Define a regular expression for the specified format
        // @todo Generalize regex to catch cases such as "foo()" and "foo" (without parentheses)	
        // Catch "foo(bar)", "foo()" and "foo"
        const std::regex pattern(R"(\s*(\w+)\((.*)\)\s*|\s*(\w+)\s*)");

        // Match the input string against the pattern
        if (std::smatch matches; std::regex_match(input, matches, pattern)) {
            // Check if there are matched groups
            if (matches.size() == 4 && (matches[1].matched || matches[3].matched)) {
                // Extract token and value from the matched groups
                token = (matches[1].matched) ? matches[1].str() : matches[3].str();
                value = matches[2].str();
                return true;
            }
        }

        // If no match is found or the match does not have the expected groups, return false
        return false;
    }
}