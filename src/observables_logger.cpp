#include "observables_logger.h"
#include "observables/observable.h"
#include "common.h"
#include "mpi.h"

#include <format>
#include <ranges>

// Constructor opens the file and writes the header
ObservablesLogger::ObservablesLogger(const std::string& filename, int bead, const std::vector<std::shared_ptr<Observable>>& observables)
    : m_bead(bead), m_observables(observables) {

    if (bead == 0) {
        m_file.open(std::format("{}/{}", Output::FOLDER_NAME, filename), std::ios::out | std::ios::app);

        if (!m_file.is_open()) {
            throw std::ios_base::failure(std::format("Failed to open {}.", Output::MAIN_FILENAME));
        }

        // Write the header
        m_file << std::format("{:^16s}", "step");

        for (const auto& observable : observables) {
            for (const auto& key : observable->quantities | std::views::keys) {
                m_file << std::vformat(" {:^16s}", std::make_format_args(key));
            }
        }

        m_file << '\n';
    }
}

// Destructor closes the file if open
ObservablesLogger::~ObservablesLogger() {
    if (m_bead == 0 && m_file.is_open()) {
        m_file.close();
    }
}

// Log observables data to the file
void ObservablesLogger::log(int step) {
    if (m_bead == 0) {
        m_file << std::format("{:^16.8e}", static_cast<double>(step));
    }

    for (const auto& observable : m_observables) {
        // The inner loop is necessary because some observable classes can calculate
        // more than one observable (e.g., "energy" calculates both the kinetic and potential energies).
        for (const double& val : observable->quantities | std::views::values) {
            double quantity_value = 0.0;
            double local_quantity_value = val;

            // Sum the results from all processes (beads)
            MPI_Allreduce(&local_quantity_value, &quantity_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            if (m_bead == 0) {
                m_file << std::format(" {:^16.8e}", quantity_value);
            }
        }
    }

    if (m_bead == 0) {
        m_file << '\n';
    }
}