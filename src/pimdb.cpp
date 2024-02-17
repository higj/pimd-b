#include "params.h"
#include "simulation.h"
#include <cstring>

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int rank, size;

	double sim_exec_time_start = 0.0;
	double sim_exec_time_end = 0.0;
	bool info_flag = false;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	std::string config_filename = "config.ini";

	try {
		for (int i = 1; i < argc; ++i) {
			if (std::strcmp(argv[i], "--dim") == 0) {
				printInfo(std::format("Program was compiled for {}-dimensional systems", NDIM), info_flag, rank);
			}
			else if (std::strcmp(argv[i], "-in") == 0) {
				// Check if there is another argument after "-in"
				if (i + 1 < argc) {
					config_filename = argv[i + 1];
					// Increment i to skip the next argument as it has been consumed as the filename
					++i;
				}
				else {
					throw std::invalid_argument("-in option requires a filename argument");
				}
			}
		}

		if (!info_flag) {
			printMsg(LOGO, rank);

			printStatus("Initializing the simulation parameters", rank);

			Params params(config_filename);

			unsigned int initial_seed = static_cast<unsigned int>(time(nullptr));

			Simulation sim(rank, size, params, initial_seed);

			printStatus("Running the simulation", rank);

			MPI_Barrier(MPI_COMM_WORLD);
			sim_exec_time_start = MPI_Wtime();

			sim.run();

			MPI_Barrier(MPI_COMM_WORLD);
			sim_exec_time_end = MPI_Wtime();

			printStatus(std::format("Simulation finished running successfully (Runtime = {:.3} sec)", sim_exec_time_end - sim_exec_time_start), rank);
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