#include "observables_logger.h"
#include "observables/observable.h"
#include "output_paths.h"
#include "mpi.h"

#include <format>
#include <ranges>

// Constructor opens the file and writes the header
ObservablesLogger::ObservablesLogger(
    const std::string& filename,
    int this_bead,
    long frequency,
    const std::vector<std::shared_ptr<Observable>>& observables
) :
    m_this_bead(this_bead),
    m_frequency(frequency),
    m_observables(observables)
{
    if (this_bead == 0)
    {
        m_file.open(std::format("{}/{}", Output::FOLDER_NAME, filename), std::ios::out | std::ios::app);

        if (!m_file.is_open())
        {
            throw std::ios_base::failure(std::format("Failed to open {}.", Output::MAIN_FILENAME));
        }

        // Write the header
        m_file << std::format("{:^16s}", "step");

        for (const auto& observable : observables)
        {
            for (const auto& key : observable->quantities | std::views::keys)
            {
                m_file << std::vformat(" {:^16s}", std::make_format_args(key));
            }
        }

        m_file << '\n';
    }
}

// Destructor closes the file if open
ObservablesLogger::~ObservablesLogger()
{
    if (m_this_bead == 0 && m_file.is_open())
    {
        m_file.close();
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
        if (m_this_bead == 0) m_file << '\n';
    //}
}

void ObservablesLogger::writeTimeStep(long step) {
    if (m_this_bead == 0) {
        m_file << std::format("{:^16.8e}", static_cast<double>(step));
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
    for (const auto& observable : m_observables) {
        for (size_t i = 0; i < observable->quantities.size(); ++i) {
            double quantity_value = global_values[idx++];

            // The simulation should be interrupted if an observable produced an invalid value
            if (!std::isfinite(quantity_value)) {
                throw std::overflow_error(
                    std::format("Invalid value of observable {}", observable->name())
                );
            }

            if (m_this_bead == 0) {
                m_file << std::format(" {:^16.8e}", quantity_value);
            }
        }
    }
}
