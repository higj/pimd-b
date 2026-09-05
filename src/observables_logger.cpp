#include "observables_logger.h"

#include "observables/observable.h"
#include "output_paths.h"
#include "mpi.h"

#include <format>
#include <map>
#include <ranges>

// Constructor opens the file and writes the header
ObservablesLogger::ObservablesLogger(
    int this_bead,
    long frequency,
    const std::vector<std::shared_ptr<Observable>>& observables
) :
    m_this_bead(this_bead),
    m_frequency(frequency),
    m_observables(observables)
{
    // Wire the cache into all observables so they can use it to store intermediate results
    for (const auto& observable : m_observables) {
        observable->setCache(&m_cache);
    }

    // Allocate MPI buffers once - observable list never changes between runs
    for (const auto& obs : m_observables) {
        m_total_quantities += static_cast<int>(obs->quantities.size());
    }
    m_local_values.resize(m_total_quantities);
    m_global_values.resize(m_total_quantities);
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

void ObservablesLogger::openFile(const std::filesystem::path& filename) {
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
        m_obs_output_file_indices.clear();
    }

    openFileAndWriteHeader(filename);
}

// Calculate and log observables data to the file
void ObservablesLogger::log(const long step)
{
    // Wipe the cache at the beginning of each logging step to ensure that observables recalculate their values
    m_cache.invalidate();

    // Calculate all observables
    for (const auto& observable : m_observables) {
        observable->calculate();
    }

    // Write the current step and the calculated observables to the output file(s)
    writeTimeStep(step);
    writeObservables();
    if (m_this_bead == 0)
    {
        for (auto& output_file : m_output_files)
        {
            output_file.stream << '\n';
        }
    }
}

/*
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
*/

void ObservablesLogger::writeTimeStep(long step) {
    if (m_this_bead == 0) {
        for (auto& output_file : m_output_files)
        {
            output_file.stream << std::format("{:^16.8e}", static_cast<double>(step));
        }
    }
}

/*
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
*/

void ObservablesLogger::writeObservables() {
    if (m_total_quantities == 0) return;

    // Fill pre-allocated send buffer
    int fill_idx = 0;
    for (const auto& observable : m_observables) {
        for (const double& val : observable->quantities | std::views::values) {
            m_local_values[fill_idx++] = val;
        }
    }

    // Sum contributions from all beads into the pre-allocated receive buffer
    MPI_Allreduce(
        m_local_values.data(),
        m_global_values.data(),
        m_total_quantities,
        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD
    );

    int idx = 0;
    for (std::size_t obs_idx = 0; obs_idx < m_observables.size(); ++obs_idx) {
        const auto& observable = m_observables[obs_idx];
        for (std::size_t q = 0; q < observable->quantities.size(); ++q) {
            const double quantity_value = m_global_values[idx++];

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
    // Here we just set up the output files and write the headers. Only rank 0 actually opens the files and writes to them
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
                m_output_files.push_back(
                    OutputFile{ .filename = output_path, .stream = {}, .observables = {} }
                );
            }

            const auto file_index = it->second;
            m_output_files[file_index].observables.push_back(observable);
            m_obs_output_file_indices.push_back(file_index);
        }

        for (auto& [out_filename, out_stream, out_observables] : m_output_files) {
            std::filesystem::create_directories(out_filename.parent_path());
            out_stream.open(out_filename, std::ios::out | std::ios::app);

            if (!out_stream.is_open()) {
                throw std::ios_base::failure(
                    std::format("Failed to open {}.", out_filename.string())
                );
            }

            out_stream << std::format("{:^16s}", "step");
            for (const auto& observable : out_observables) {
                for (const auto& key : observable->quantities | std::views::keys) {
                    out_stream << std::vformat(" {:^16s}", std::make_format_args(key));
                }
            }
            out_stream << '\n';
        }
    }
}