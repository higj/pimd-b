#pragma once

#include "observables/observable_cache.h"

#include <memory>
#include <vector>
#include <fstream>
#include <filesystem>

class Observable;

class ObservablesLogger {
public:
    ObservablesLogger(
        const std::filesystem::path& filename,
        int this_bead,
        long frequency,
        const std::vector<std::shared_ptr<Observable>>& observables
    );
    ~ObservablesLogger();

    bool isLoggingStep(const long step) const { return step % m_frequency == 0; }

    void log(long step);

    /**
     * @brief Close the current output file and reopen a new one with the specified filename.
     * This is useful for RPMD mode where each run needs a separate output file.
     * @param new_filename The name of the new file to open for logging observables.
     */
    void reopenFile(const std::filesystem::path& new_filename);

private:
    struct OutputFile {
        std::filesystem::path filename;
        std::ofstream stream;
        std::vector<std::shared_ptr<Observable>> observables;
    };

    int m_this_bead;
    long m_frequency;
    std::vector<std::shared_ptr<Observable>> m_observables;

    std::filesystem::path m_main_output_filename;
    std::vector<OutputFile> m_output_files;
    /// Maps each observable (by index into m_observables) to its OutputFile index in m_output_files.
    std::vector<std::size_t> m_obs_output_file_indices;

    ObservableCache m_cache;

    std::vector<double> m_local_values;
    std::vector<double> m_global_values;
    int m_total_quantities = 0;

    void writeTimeStep(long step);
    void writeObservables();
    void openFileAndWriteHeader(const std::filesystem::path& filename);
};
