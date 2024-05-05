#include "params.h"
#include "simulation.h"
#include <cstring>
#include "mpi.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::string config_filename = "config.ini";

    try {
        bool info_flag = false;

        // Check for command line arguments.
        // If the "--dim" flag is present, print the dimension of the system and exit.
        // If the "-in" flag is present, use the next argument as the configuration filename. Otherwise, use the default filename.
        for (int i = 1; i < argc; ++i) {
            if (std::strcmp(argv[i], "--dim") == 0) {
                printInfo(std::format("Program was compiled for {}-dimensional systems", NDIM), info_flag, rank);
            } else if (std::strcmp(argv[i], "--bosonic_alg") == 0) {
                if constexpr (OLD_BOSONIC_ALGORITHM)
                    printInfo("Program was compiled with original bosonic algorithm.", info_flag, rank);
                else
                    printInfo("Program was compiled with quadratic bosonic algorithm.", info_flag, rank);
            } else if (std::strcmp(argv[i], "-in") == 0) {
                // Check if there is another argument after "-in"
                if (i + 1 < argc) {
                    config_filename = argv[i + 1];
                    // Increment i to skip the next argument as it has been consumed as the filename
                    ++i;
                } else {
                    throw std::invalid_argument("-in option requires a filename argument");
                }
            } else {
                throw std::invalid_argument(std::format("Unknown option: {}", argv[i]));
            }
        }

        // If we got to this point, and no info has been requested then initiate the simulation
        if (!info_flag) {
            printMsg(LOGO, rank);
            printStatus("Initializing the simulation parameters", rank);

            // Load the simulation parameters from the configuration file
            Params params(config_filename);

            // Initialize the random number generator seed using the current time
            const unsigned int initial_seed = static_cast<unsigned int>(time(nullptr));

            Simulation sim(rank, size, params, initial_seed);

            printStatus("Running the simulation", rank);

            MPI_Barrier(MPI_COMM_WORLD);
            const double sim_exec_time_start = MPI_Wtime();

            sim.run();

            MPI_Barrier(MPI_COMM_WORLD);
            const double sim_exec_time_end = MPI_Wtime();

            printStatus(std::format("Simulation finished running successfully (Runtime = {:.3} sec)",
                                    sim_exec_time_end - sim_exec_time_start), rank);
        }
    }
    catch (const std::invalid_argument& ex) {
        printError(ex.what(), rank, ErrorMessage::INVALID_ARG_ERR);
    }
    catch (const std::overflow_error& ex) {
        printError(ex.what(), rank, ErrorMessage::OVERFLOW_ERR);
    }
    catch (const std::exception& ex) {
        printError(ex.what(), rank, ErrorMessage::GENERAL_ERR);
    }

    MPI_Finalize();

    return 0;
}
