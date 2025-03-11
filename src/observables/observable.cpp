#include "observables.h"
#include "simulation.h"
#include "units.h"
#include <ranges>
#include "mpi.h"

/**
 * @brief Generic observable class constructor
 */
Observable::Observable(const Simulation& _sim, int _freq, const std::string& _out_unit) :
    sim(_sim), freq(_freq), out_unit(_out_unit) {
}

/**
 * Initializes observables with the given labels.
 *
 * @param labels Labels of the quantities to be calculated
 */
void Observable::initialize(const std::vector<std::string>& labels) {
    for (const std::string& label : labels) {
        quantities.insert({ label, 0.0 });
    }
}

/**
 * @brief Resets the values of the observable to zero.
 * Useful for clearing results from the previous moleclar dynamics step.
 */
void Observable::resetValues() {
    for (auto it = quantities.begin(); it != quantities.end(); ++it) {
        it.value() = 0.0;
    }
}

/**
 * @brief Creates an observable object based on the observable type.
 *
 * @param observable_type Type of observable to create
 * @param _sim Simulation object
 * @param _freq Frequency at which the observable is calculated
 * @param _out_unit Output unit of the observable
 * @return std::unique_ptr<Observable> Pointer to the created observable object
 * @throw std::invalid_argument If the observable type is unknown
 */
std::unique_ptr<Observable> ObservableFactory::createQuantity(const std::string& observable_type,
    const Simulation& _sim, int _freq,
    const std::string& _out_unit) {
    if (observable_type == "energy") {
        return std::make_unique<EnergyObservable>(_sim, _freq, _out_unit);
    } else if (observable_type == "classical") {
        return std::make_unique<ClassicalObservable>(_sim, _freq, _out_unit);
    } else if (observable_type == "bosonic") {
        return std::make_unique<BosonicObservable>(_sim, _freq, _out_unit);
    } else if (observable_type == "gsf") {
        return std::make_unique<GSFActionObservable>(_sim, _freq, _out_unit);
    } else {
        throw std::invalid_argument("Unknown observable type.");
    }
}

// Constructor opens the file and writes the header
ObservablesLogger::ObservablesLogger(const std::string& filename, int _bead, const std::vector<std::unique_ptr<Observable>>& _observables)
    : bead(_bead), observables(_observables) {

    if (bead == 0) {
        file.open(std::format("{}/{}", Output::FOLDER_NAME, filename), std::ios::out | std::ios::app);

        if (!file.is_open()) {
            throw std::ios_base::failure(std::format("Failed to open {}.", Output::MAIN_FILENAME));
        }

        // Write the header
        file << std::format("{:^16s}", "step");

        for (const auto& observable : observables) {
            for (const auto& key : observable->quantities | std::views::keys) {
                file << std::vformat(" {:^16s}", std::make_format_args(key));
            }
        }

        file << '\n';
    }    
}

// Destructor closes the file if open
ObservablesLogger::~ObservablesLogger() {
    if (bead == 0 && file.is_open()) {
        file.close();
    }
}

// Log observables data to the file
void ObservablesLogger::log(int step) {
    if (bead == 0) {
        file << std::format("{:^16.8e}", static_cast<double>(step));
    }

    for (const auto& observable : observables) {
        // The inner loop is necessary because some observable classes can calculate
        // more than one observable (e.g., "energy" calculates both the kinetic and potential energies).
        for (const double& val : observable->quantities | std::views::values) {
            double quantity_value = 0.0;
            double local_quantity_value = val;

            // Sum the results from all processes (beads)
            MPI_Allreduce(&local_quantity_value, &quantity_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            if (bead == 0) {
                file << std::format(" {:^16.8e}", quantity_value);
            }
        }
    }

    if (bead == 0) {
        file << '\n';
    }
}