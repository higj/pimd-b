#pragma once

#include <memory>
#include <vector>
#include <fstream>

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

    bool isLoggingStep(long step) const { return step % m_frequency == 0; }

    void log(long step);

    /**
     * @brief Close the current output file and reopen a new one with the specified filename.
     * This is useful for RPMD mode where each run needs a separate output file.
     * @param new_filename The name of the new file to open for logging observables.
     */
    void reopenFile(const std::filesystem::path& new_filename);

private:
    std::ofstream m_file;
    int m_this_bead;
    long m_frequency;
    std::vector<std::shared_ptr<Observable>> m_observables;

    void writeTimeStep(long step);
    void writeObservables();
    void openFileAndWriteHeader(const std::filesystem::path& filename);
};