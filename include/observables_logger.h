#pragma once

#include "observables/observable_cache.h"

#include <memory>
#include <vector>
#include <fstream>
#include <filesystem>

class Observable;

class ObservablesLogger {
public:
    /**
     * @brief Wires cache into observables, allocates MPI buffers.
     * No file is opened; all ranks must call this.
     *
     * @param this_bead Current bead index (0-based) for the calling MPI rank.
     * @param frequency Frequency of logging steps (every N steps).
     * @param observables Vector of shared pointers to the observables to be logged.
     */
    ObservablesLogger(
        int this_bead,
        long frequency,
        const std::vector<std::shared_ptr<Observable>>& observables
    );

    ~ObservablesLogger();

    /**
     * @brief Opens (or re-opens) the output file and writes the column headers.
     * Call once before the first step of a PIMD run, or at the start of each
     * NVE run of an RPMD simulation. All ranks must call this, but only
     * rank 0 actually opens the file. The MPI buffer sizes are re-confirmed here,
     * but the resize is a no-op after the first call since the observable list
     * never changes.
     *
     * @param filename The name of the file to open for logging observables. This is relative to the active simulation output directory.
     */
    void openFile(const std::filesystem::path& filename);

    /**
     * @brief Check if the given step is a logging step.
     *
     * @param step The current simulation step.
     * @return True if the step is a logging step, false otherwise.
     */
    [[nodiscard]] bool isLoggingStep(const long step) const { return step % m_frequency == 0; }

    /**
     * @brief Log the current values of all observables to the output file(s).
     *
     * @param step The current simulation step. This is used to write the step number in the output file.
     */
    void log(long step);

    /**
     * @brief Close the current output file and reopen a new one with the specified filename.
     * This is useful for RPMD mode where each run needs a separate output file.
     * @param new_filename The name of the new file to open for logging observables.
     */
    //void reopenFile(const std::filesystem::path& new_filename);

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
