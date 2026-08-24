#include <cstring>
#include <iostream>
#include <string_view>

#include "mpi.h"
#include "simulation.h"
#include "app_config.h"
#include "error_messages.h"
#include "params.h"

namespace {
    // Prints one of two compile-time-fixed messages depending on `flag`, and marks
    // that the user asked for info (as opposed to requesting a simulation run).
    void printCompileTimeFlag(bool flag, const char* on_msg, const char* off_msg, bool& info_flag, int rank) {
        printMsg(flag ? on_msg : off_msg, rank);
        info_flag = true;
    }

    void parseArguments(const int arg_num, char** arg_arr, std::string& conf_filename, bool& info_flag, const int rank) {
        // Check for command line arguments.
        // If the "--dim" flag is present, print the dimension of the system and exit.
        // If the "--bosonic_alg" flag is present, print the bosonic algorithm used and exit.
        // If the "--convention" flag is present, print the convention used and exit.
        // If the "-in" flag is present, use the next argument as the configuration filename. Otherwise, use the default filename.
        for (int i = 1; i < arg_num; ++i) {
            const std::string_view arg{ arg_arr[i] };

            if (arg == "--dim") {
                printMsg(std::format("Program was compiled for {}-dimensional systems", NDIM), rank);
                info_flag = true;
            } else if (arg == "--bosonic_alg") {
                printCompileTimeFlag(
                    FACTORIAL_BOSONIC_ALGORITHM,
                    "Program was compiled with naive bosonic algorithm.",
                    "Program was compiled with quadratic bosonic algorithm.",
                    info_flag, rank);
            } else if (arg == "--convention") {
                printCompileTimeFlag(
                    TAU_CONVENTION,
                    "Program was compiled with tau convention.",
                    "Program was compiled with beta convention.",
                    info_flag, rank);
            } else if (arg == "-in") {
                // Check if there is another argument after "-in"
                if (i + 1 < arg_num) {
                    conf_filename = arg_arr[i + 1];
                    // Increment i to skip the next argument as it has been consumed as the filename
                    ++i;
                } else {
                    throw std::invalid_argument("-in option requires a filename argument");
                }
            } else {
                throw std::invalid_argument(
                    std::format("Unrecognized command line argument: '{}'", arg));
            }
        }
    }

    void throwError(
        const std::exception& ex,
        int& error_flag,
        int rank,
        const std::string& error_type
    ) {
        printError(ex.what(), rank, error_type);
        error_flag = 1;
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int error_flag = 0;  // 0 = no error, 1 = error occurred

    try {
        // Flag to check if the user requested information about the program as opposed to running the simulation
        bool display_info = false;

        std::string config_filename = AppConfig::DEFAULT_CONFIG_FILENAME;
        parseArguments(argc, argv, config_filename, display_info, rank);

        // If we got to this point, and no info has been requested then initiate the simulation
        if (!display_info) {
            printMsg(AppConfig::LOGO, rank);

            // Parse the simulation parameters from the configuration (input) file and run the simulation
            const auto config = Params(config_filename, rank).load();
            Simulation sim(rank, size, config);
            sim.run();
        }
    } catch (const std::invalid_argument& ex) {
        throwError(ex, error_flag, rank, ErrorMessage::INVALID_ARG_ERR);
    } catch (const std::overflow_error& ex) {
        throwError(ex, error_flag, rank, ErrorMessage::OVERFLOW_ERR);
    } catch (const std::exception& ex) {
        throwError(ex, error_flag, rank, ErrorMessage::GENERAL_ERR);
    }

    // Flush output buffers to ensure messages are displayed before MPI operations
    std::cout.flush();
    std::cerr.flush();
    fflush(stdout);
    fflush(stderr);

    // Broadcast error status to ensure all ranks know about the error
    int global_error_flag = 0;
    MPI_Allreduce(&error_flag, &global_error_flag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    if (global_error_flag != 0) {
        // Error occurred on one or more ranks - terminate all ranks cleanly
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}