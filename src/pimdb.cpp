#include <cstring>
#include "mpi.h"
#include "params.h"
#include "simulation.h"

void parseArguments(int arg_num, char** arg_arr, std::string& conf_filename, bool& info_flag, const int rank) {
    // Check for command line arguments.
    // If the "--dim" flag is present, print the dimension of the system and exit.
    // If the "-in" flag is present, use the next argument as the configuration filename. Otherwise, use the default filename.
    for (int i = 1; i < arg_num; ++i) {
        if (std::strcmp(arg_arr[i], "--dim") == 0) {
            printInfo(std::format("Program was compiled for {}-dimensional systems", NDIM), info_flag, rank);
        } else if (std::strcmp(arg_arr[i], "--bosonic_alg") == 0) {
            if constexpr (OLD_BOSONIC_ALGORITHM)
                printInfo("Program was compiled with naive bosonic algorithm.", info_flag, rank);
            else
                printInfo("Program was compiled with quadratic bosonic algorithm.", info_flag, rank);
        } else if (std::strcmp(arg_arr[i], "-in") == 0) {
            // Check if there is another argument after "-in"
            if (i + 1 < arg_num) {
                conf_filename = arg_arr[i + 1];
                // Increment i to skip the next argument as it has been consumed as the filename
                ++i;
            } else {
                throw std::invalid_argument("-in option requires a filename argument");
            }
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::string config_filename = "config.ini";

    try {
        // Flag to check if the user requested information about the program as opposed to running the simulation
        bool display_info = false;

        parseArguments(argc, argv, config_filename, display_info, rank);

        // If we got to this point, and no info has been requested then initiate the simulation
        if (!display_info) {
            printMsg(LOGO, rank);

            // Load the simulation parameters from the configuration file
            Params params(config_filename, rank);
            // Initialize the random number generator seed based on the current time
            Simulation sim(rank, size, params, static_cast<unsigned int>(time(nullptr)));
            sim.run();
        }
    } catch (const std::invalid_argument& ex) {
        printError(ex.what(), rank, ErrorMessage::INVALID_ARG_ERR);
    } catch (const std::overflow_error& ex) {
        printError(ex.what(), rank, ErrorMessage::OVERFLOW_ERR);
    } catch (const std::exception& ex) {
        printError(ex.what(), rank, ErrorMessage::GENERAL_ERR);
    }

    MPI_Finalize();

    return 0;
}
