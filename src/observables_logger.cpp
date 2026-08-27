#include "observables_logger.h"

#include "observables/observable.h"
#include "output_paths.h"
#include "mpi.h"

#include <format>
#include <map>
#include <ranges>

// Constructor opens the file and writes the header
ObservablesLogger::ObservablesLogger(
    const std::filesystem::path& filename,
    int this_bead,
    long frequency,
    const std::vector<std::shared_ptr<Observable>>& observables
) :
    m_this_bead(this_bead),
    m_frequency(frequency),
    m_observables(observables)
{
    openFileAndWriteHeader(filename);
}

// Destructor closes the file if open
ObservablesLogger::~ObservablesLogger()
{
    if (m_this_bead == 0)
    {
        for (auto& output_file : m_output_files)
        {
            if (output_file.stream.is_open())
            {
                output_file.stream.close();
            }
        }
    }
}

// Log observables data to the file
void ObservablesLogger::log(long step)
{
    /// TODO: Better do this check in Simulation::calculateAndLogObservables, before calculating observables 
    ///       (if we won't write them, why bother calculating?)
    //if (step % m_frequency == 0)
    //{
        writeTimeStep(step);
        writeObservables();
        if (m_this_bead == 0)
        {
            for (auto& output_file : m_output_files)
            {
                output_file.stream << '\n';
            }
        }
    //}
}

void ObservablesLogger::reopenFile(const std::filesystem::path& new_filename) {
    if (m_this_bead == 0)
    {
        for (auto& output_file : m_output_files)
        {
            if (output_file.stream.is_open())
            {
                output_file.stream.close();
            }
        }
        m_output_files.clear();

        openFileAndWriteHeader(new_filename);
    }
}

void ObservablesLogger::writeTimeStep(long step) {
    if (m_this_bead == 0) {
        for (auto& output_file : m_output_files)
        {
            output_file.stream << std::format("{:^16.8e}", static_cast<double>(step));
        }
    }
}

void ObservablesLogger::writeObservables() {
    std::vector<double> local_values;
    for (const auto& observable : m_observables) {
        // The inner loop is necessary because some observable classes can calculate
        // more than one observable (e.g., "energy" calculates both the kinetic and potential energies).
        for (const double& val : observable->quantities | std::views::values) {
            local_values.push_back(val);
        }
    }

    if (local_values.empty()) return;

    std::vector<double> global_values(local_values.size(), 0.0);
    // Sum the results from all processes (beads)
    MPI_Allreduce(local_values.data(), global_values.data(), static_cast<int>(local_values.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    int idx = 0;
    for (std::size_t obs_idx = 0; obs_idx < m_observables.size(); ++obs_idx) {
        const auto& observable = m_observables[obs_idx];
        for (std::size_t i = 0; i < observable->quantities.size(); ++i) {
            const double quantity_value = global_values[idx++];

            // The simulation should be interrupted if an observable produced an invalid value
            if (!std::isfinite(quantity_value)) {
                throw std::overflow_error(
                    std::format("Invalid value of observable {}", observable->name())
                );
            }

            if (m_this_bead == 0) {
                const std::size_t file_idx = m_obs_output_file_indices[obs_idx];
                m_output_files[file_idx].stream
                    << std::format(" {:^16.8e}", quantity_value);
            }
        }
    }
}

void ObservablesLogger::openFileAndWriteHeader(const std::filesystem::path& filename) {
    if (m_this_bead == 0) {
        std::map<std::filesystem::path, std::size_t> file_indices;

        m_output_files.clear();
        m_obs_output_file_indices.clear();

        m_main_output_filename = filename;

        for (const auto& observable : m_observables) {
            const auto output_path = observable->usesCustomFile()
                ? filename.parent_path() / observable->outputFilename()
                : filename;

            const auto [it, inserted] = file_indices.emplace(
                output_path, m_output_files.size()
            );

            if (inserted) {
                m_output_files.push_back(OutputFile{ .filename = output_path });
            }

            const auto file_index = it->second;
            m_output_files[file_index].observables.push_back(observable);
            m_obs_output_file_indices.push_back(file_index);
        }

        // Custom names are relative to the active output directory; the default
        // continues to use Simulation's supplied filename (including RPMD paths).
        for (auto& output_file : m_output_files)
        {
            std::filesystem::create_directories(output_file.filename.parent_path());
            output_file.stream.open(output_file.filename, std::ios::out | std::ios::app);

            if (!output_file.stream.is_open()) {
                throw std::ios_base::failure(
                    std::format("Failed to open {}.", output_file.filename.string())
                );
            }

            output_file.stream << std::format("{:^16s}", "step");
            for (const auto& observable : output_file.observables) {
                for (const auto& key : observable->quantities | std::views::keys) {
                    output_file.stream << std::vformat(" {:^16s}", std::make_format_args(key));
                }
            }

            output_file.stream << '\n';
        }
    }
}
