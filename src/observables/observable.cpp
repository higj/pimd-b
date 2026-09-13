#include "observables.h"
#include "output_paths.h"

#include <filesystem>

Observable::Observable(
    const std::string& name, 
    int out_freq, 
    const std::string& out_unit,
    const std::string& out_filename
) : m_name(name), m_out_freq(out_freq), m_out_unit(out_unit), m_out_filename(out_filename)
{
    initializeFolder(Output::FOLDER_NAME);
}

void Observable::initializeFolder(const std::string_view folder_name)
{
    std::filesystem::create_directory(folder_name);
}

void Observable::initializeLabel(const std::string& label)
{
    quantities.insert({label, 0.0});
}

void Observable::initialize(const std::vector<std::string>& labels)
{
    for (const std::string& label : labels)
    {
        initializeLabel(label);
    }
}

void Observable::resetValues()
{
    for (auto it = quantities.begin(); it != quantities.end(); ++it)
    {
        it.value() = 0.0;
    }
}
