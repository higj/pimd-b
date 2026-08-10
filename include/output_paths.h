#pragma once

#include <string>
#include <filesystem>

namespace Output {
    inline const std::string FOLDER_NAME = "output";
    inline const std::string MAIN_FILENAME = "simulation.out";

    inline std::filesystem::path getPimdFilename()
    {
        return std::filesystem::path(FOLDER_NAME) / MAIN_FILENAME;
    }

    inline std::filesystem::path getRpmdFolder(const int run) {
        return std::filesystem::path(FOLDER_NAME) / ("rpmd_" + std::to_string(run));
    }

    /**
     * @brief Generate the file path for the RPMD output file based on the run number.
     * @param run The sub-simulation run number (0-based index).
     * @return File path in the format "rpmd_<run>/simulation.out".
     */
    inline std::filesystem::path getRpmdFilename(const int run) {
        // Construct the relative path
        return getRpmdFolder(run) / MAIN_FILENAME;

        // Create rpmd_<run>.out file in the current working directory
        // Disabled, because we want to create a directory for each run instead of just a file in the current directory
        // return std::string("rpmd_") + std::to_string(run) + ".out";

        // Return the full file path as a string
        //return filepath.string();
    }
}
