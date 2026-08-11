#pragma once

#include <string>
#include <filesystem>

namespace Output {
    inline const std::string FOLDER_NAME = "output";
    inline const std::string MAIN_FILENAME = "simulation.out";

    /**
     * @brief Generate the file path for the PIMD output file
     * @return File path in the format "output/simulation.out".
     */
    inline std::filesystem::path getPimdFilename()
    {
        return std::filesystem::path(FOLDER_NAME) / MAIN_FILENAME;
    }

    /**
     * @brief Generate the file path for the RPMD output folder based on the run number.
     * @param run The sub-simulation run number (0-based index).
     * @return File path in the format "rpmd_<run>".
     */
    inline std::filesystem::path getRpmdFolder(const int run) {
        return std::filesystem::path(FOLDER_NAME) / ("rpmd_" + std::to_string(run));
    }

    /**
     * @brief Generate the file path for the RPMD output file based on the run number.
     * @param run The sub-simulation run number (0-based index).
     * @return File path in the format "rpmd_<run>/simulation.out".
     */
    inline std::filesystem::path getRpmdFilename(const int run) {
        return getRpmdFolder(run) / MAIN_FILENAME;
    }
}
