#include "params.h"
#include "simulation.h"
#include <cstring>

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int rank, size;
	double sim_exec_time_start, sim_exec_time_end;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	for (int i = 1; i < argc; ++i) {
		if (std::strcmp(argv[i], "--dim") == 0) {
			std::cout << "Program was compiled for " << NDIM << "-dimensional systems." << std::endl;

			return 0;
		}
	}

	printMsg(LOGO, rank);

	try {
		printStatus("Initializing the simulation parameters", rank);

		Params params;

		unsigned int initial_seed = static_cast<unsigned int>(time(nullptr));

		Simulation sim(rank, size, params, initial_seed);

		printStatus("Running the simulation", rank);

		MPI_Barrier(MPI_COMM_WORLD);
		sim_exec_time_start = MPI_Wtime();

		sim.run();

		MPI_Barrier(MPI_COMM_WORLD);
		sim_exec_time_end = MPI_Wtime();
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

	printStatus(std::format("Simulation finished running successfully (Runtime = {:.3} sec)", sim_exec_time_end - sim_exec_time_start), rank);

	return 0;
}