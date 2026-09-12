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
    if (m_this_bead != 0) return;

#ifdef USE_HDF5
    for (auto& h5_f : m_h5_obs_files) {
        for (const hid_t ds : h5_f.col_datasets) {
            if (ds != H5I_INVALID_HID) {
                H5Dclose(ds);
            }
        }

        if (h5_f.step_ds != H5I_INVALID_HID) {
            H5Dclose(h5_f.step_ds);
        }
        if (h5_f.file_id != H5I_INVALID_HID) {
            H5Fclose(h5_f.file_id);
        }
    }
#else
    for (auto& output_file : m_output_files)
    {
        if (output_file.stream.is_open())
        {
            output_file.stream.close();
        }
    }
#endif
}

void ObservablesLogger::openFile(const std::filesystem::path& filename) {
    if (m_this_bead == 0)
    {
#ifdef USE_HDF5
        // Close any previously open HDF5 files
        for (auto& h5_f : m_h5_obs_files) {
            for (const hid_t ds : h5_f.col_datasets)
                if (ds != H5I_INVALID_HID) H5Dclose(ds);
            if (h5_f.step_ds != H5I_INVALID_HID) H5Dclose(h5_f.step_ds);
            if (h5_f.file_id != H5I_INVALID_HID) H5Fclose(h5_f.file_id);
        }
        m_h5_obs_files.clear();
        m_h5_obs_file_indices.clear();
        m_h5_obs_col_offsets.clear();
#else
        for (auto& output_file : m_output_files)
        {
            if (output_file.stream.is_open())
            {
                output_file.stream.close();
            }
        }

        m_output_files.clear();
        m_obs_output_file_indices.clear();
#endif
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
#ifdef USE_HDF5
        // Advance the row count for each HDF5 file so that the next logging step writes to the next row
        for (auto& h5_f : m_h5_obs_files) {
            ++h5_f.row_count;
        }
#else
        for (auto& output_file : m_output_files)
        {
            output_file.stream << '\n';
        }
#endif
    }
}

void ObservablesLogger::writeTimeStep(long step)
{
    if (m_this_bead != 0) return;

#ifdef USE_HDF5
    for (const auto& h5_f : m_h5_obs_files) {
        H5Utils::append_int(h5_f.step_ds, h5_f.row_count, static_cast<int>(step));
    }
#else
    for (auto& output_file : m_output_files) {
        output_file.stream << std::format("{:^16.8e}", static_cast<double>(step));
    }
#endif
}

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

    if (m_this_bead != 0) return; // Only rank 0 writes to the output files

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

#ifdef USE_HDF5
            const std::size_t file_idx = m_h5_obs_file_indices[obs_idx];
            auto& h5_f = m_h5_obs_files[file_idx];
            const std::size_t col = m_h5_obs_col_offsets[obs_idx] + q;
            H5Utils::append_double(h5_f.col_datasets[col], h5_f.row_count, quantity_value);
#else
            const std::size_t file_idx = m_obs_output_file_indices[obs_idx];
            m_output_files[file_idx].stream
                << std::format(" {:^16.8e}", quantity_value);
#endif
        }
    }
}

void ObservablesLogger::openFileAndWriteHeader(const std::filesystem::path& filename) {
    // Here we just set up the output files and write the headers. Only rank 0 actually opens the files and writes to them
    if (m_this_bead != 0) return;

#ifdef USE_HDF5
    // Derive the .h5 filename from the .out filename (e.g. simulation.out -> simulation.h5)
    auto h5_of = [](const std::filesystem::path& txt_path) {
        const std::string stem = txt_path.stem().string();
        return txt_path.parent_path() / (stem + ".h5");
        };

    std::map<std::filesystem::path, std::size_t> h5_file_map;
    m_main_output_filename = filename;

    for (const auto& obs : m_observables)
    {
        const auto txt_path = obs->usesCustomFile()
                                  ? filename.parent_path() / obs->outputFilename()
                                  : filename;
        const auto h5_path = h5_of(txt_path);

        auto [it, inserted] = h5_file_map.emplace(h5_path, m_h5_obs_files.size());
        if (inserted) {
            std::filesystem::create_directories(h5_path.parent_path());
            H5ObsFile h5_f;
            h5_f.file_id = H5Fcreate(h5_path.string().c_str(),
                H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT
            );

            if (h5_f.file_id < 0) {
                throw std::ios_base::failure("Failed to create HDF5 file: " + h5_path.string());
            }

            h5_f.step_ds = H5Utils::make_1d(h5_f.file_id, "step", H5T_NATIVE_INT);
            m_h5_obs_files.push_back(std::move(h5_f));
        }

        const std::size_t file_idx = it->second;
        auto& h5_f = m_h5_obs_files[file_idx];

        m_h5_obs_file_indices.push_back(file_idx);
        m_h5_obs_col_offsets.push_back(h5_f.col_datasets.size());

        // One dataset per quantity, named by the observable's quantity key
        for (const auto& key : obs->quantities | std::views::keys) {
            h5_f.col_datasets.push_back(
                H5Utils::make_1d(h5_f.file_id, key.c_str(), H5T_NATIVE_DOUBLE)
            );
        }
    }
#else
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
#endif
}
