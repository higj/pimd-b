// Read an INI file into easy-to-access name/value pairs.

// SPDX-License-Identifier: BSD-3-Clause

// Copyright (C) 2009-2020, Ben Hoyt

// inih and INIReader are released under the New BSD license (see LICENSE.txt).
// Go to the project home page for more info:
//
// https://github.com/benhoyt/inih

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include "ini.h"
#include "inireader.h"

INIReader::INIReader(const std::string& filename)
{
    _error = ini_parse(filename.c_str(), ValueHandler, this);
}

INIReader::INIReader(const char* buffer, size_t buffer_size)
{
    const std::string content(buffer, buffer_size);
    _error = ini_parse_string(content.c_str(), ValueHandler, this);
}

int INIReader::ParseError() const
{
    return _error;
}

std::string INIReader::Get(const std::string& section, const std::string& name, const std::string& default_value) const
{
    std::string key = MakeKey(section, name);
    // Use _values.find() here instead of _values.at() to support pre C++11 compilers
    return _values.contains(key) ? _values.find(key)->second : default_value;
}

std::string INIReader::GetString(const std::string& section, const std::string& name,
                                 const std::string& default_value) const
{
    const std::string str = Get(section, name, "");
    return str.empty() ? default_value : str;
}

long INIReader::GetLong(const std::string& section, const std::string& name, long default_value) const
{
    const std::string valstr = Get(section, name, "");
    const char* value = valstr.c_str();
    char* end;
    // This parses "1234" (decimal) and also "0x4D2" (hex)
    long n = strtol(value, &end, 0);
    return end > value ? n : default_value;
}

unsigned long INIReader::GetUnsigned(const std::string& section, const std::string& name, unsigned long default_value) const {
    const std::string valstr = Get(section, name, "");
    const char* value = valstr.c_str();
    char* end;
    // This parses "1234" (decimal) and also "0x4D2" (hex)
    unsigned long n = strtoul(value, &end, 0);
    return end > value ? n : default_value;
}

int INIReader::GetInteger(const std::string& section, const std::string& name, int default_value) const
{
    std::string valstr = Get(section, name, "");
    const char* value = valstr.c_str();
    int n = atoi(value);

    return (value && !*value) ? default_value : n;
}

double INIReader::GetReal(const std::string& section, const std::string& name, double default_value) const
{
    std::string valstr = Get(section, name, "");
    const char* value = valstr.c_str();
    char* end;
    double n = strtod(value, &end);
    return end > value ? n : default_value;
}

bool INIReader::GetBoolean(const std::string& section, const std::string& name, bool default_value) const
{
    std::string valstr = Get(section, name, "");
    // Convert to lower case to make string comparisons case-insensitive
    LowerString(valstr);
    if (valstr == "true" || valstr == "yes" || valstr == "on" || valstr == "1")
        return true;
    if (valstr == "false" || valstr == "no" || valstr == "off" || valstr == "0")
        return false;
    return default_value;
}

bool INIReader::HasSection(const std::string& section) const
{
    const std::string key = MakeKey(section, "");
    auto pos = _values.lower_bound(key);
    if (pos == _values.end())
        return false;
    // Does the key at the lower_bound pos start with "section"?
    return pos->first.starts_with(key);
}

bool INIReader::HasValue(const std::string& section, const std::string& name) const
{
    const std::string key = MakeKey(section, name);
    return _values.contains(key);
}

void INIReader::LowerString(std::string& str_to_lower)
{
    //std::ranges::transform(str_to_lower, str_to_lower.begin(),
    //    [](const unsigned char& ch) { return static_cast<unsigned char>(::tolower(ch)); });

    std::ranges::transform(str_to_lower, str_to_lower.begin(),
                           [](const unsigned char& c)
                           {
                               return static_cast<unsigned char>(std::tolower(c));
                           });
}

std::string INIReader::MakeKey(const std::string& section, const std::string& name)
{
    std::string key = section + "=" + name;
    // Convert to lower case to make section/name lookups case-insensitive
    LowerString(key);
    return key;
}

int INIReader::ValueHandler(void* user, const char* section, const char* name,
                            const char* value)
{
    if (!name) // Happens when INI_CALL_HANDLER_ON_NEW_SECTION enabled
        return 1;
    const auto reader = static_cast<INIReader*>(user);
    const std::string key = MakeKey(section, name);
    if (!reader->_values[key].empty())
        reader->_values[key] += "\n";
    reader->_values[key] += value ? value : "";
    return 1;
}
