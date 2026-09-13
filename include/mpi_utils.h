#pragma once

#include <mpi.h>

#include <algorithm>
#include <array>
#include <climits>
#include <format>
#include <stdexcept>
#include <string>
#include <string_view>

namespace MpiUtils {

    constexpr std::size_t ERROR_MESSAGE_CAPACITY = 1024;

    /**
     * Executes an operation collectively.
     *
     * Every MPI rank must call this function in the same order. If the operation
     * throws on any rank, all ranks receive the same exception before proceeding.
     */
    template <typename Operation>
    void runCollectively(
        std::string_view operation_name,
        Operation&& operation,
        MPI_Comm communicator = MPI_COMM_WORLD
    ) {
        int rank = 0;
        MPI_Comm_rank(communicator, &rank);

        std::string local_error;

        try {
            operation();
        } catch (const std::exception& error) {
            local_error = error.what();
        } catch (...) {
            local_error = "Unknown non-standard exception.";
        }

        const int failure_rank_candidate =
            local_error.empty() ? INT_MAX : rank;

        int failure_rank = INT_MAX;
        MPI_Allreduce(
            &failure_rank_candidate,
            &failure_rank,
            1,
            MPI_INT,
            MPI_MIN,
            communicator
        );

        if (failure_rank == INT_MAX) {
            return;
        }

        std::array<char, ERROR_MESSAGE_CAPACITY> error_buffer{};
        int error_length = 0;

        if (rank == failure_rank) {
            error_length = static_cast<int>(std::min(
                local_error.size(),
                ERROR_MESSAGE_CAPACITY - 1
            ));

            std::copy_n(
                local_error.data(),
                error_length,
                error_buffer.data()
            );
        }

        MPI_Bcast(
            &error_length,
            1,
            MPI_INT,
            failure_rank,
            communicator
        );

        MPI_Bcast(
            error_buffer.data(),
            error_length,
            MPI_CHAR,
            failure_rank,
            communicator
        );

        throw std::runtime_error(std::format(
            "{0} failed on MPI rank {1}\n{3:-^40}\n{2}\n{3:-^40}",
            operation_name,
            failure_rank,
            std::string_view(error_buffer.data(), error_length),
            std::string()
        ));
    }

}