#include <iostream>

#include "common.h"

void printMsg(const std::string& msg, int this_rank, int out_rank) {
    if (this_rank == out_rank) {
        std::cout << msg << '\n';
    }
}

void printInfo(const std::string& info, bool& info_flag, int this_rank, int out_rank) {
    printMsg(info, this_rank, out_rank);
    info_flag = true;
}

void printStatus(const std::string& status, int this_rank, int out_rank) {
    printMsg("[*] " + status, this_rank, out_rank);
}

void printError(const std::string& msg, int this_rank, const std::string& err_type, int out_rank) {
    printMsg("[X] " + err_type + ": " + msg, this_rank, out_rank);
}

void printProgress(int this_step, int total_steps, int this_rank, int out_rank) {
    if (this_rank == out_rank) {
        double percentage = static_cast<double>(this_step) / (total_steps - 1);
        int val = static_cast<int>(percentage * 100);
        int lpad = static_cast<int>(percentage * PBWIDTH);
        int rpad = PBWIDTH - lpad;

        printf("\r[%.*s%*s] %3d%%", lpad, PBSTR.data(), rpad, "", val);
        fflush(stdout);

        if (this_step == total_steps - 1) {
            printf("\n");
        }
    }
}

void applyMinimumImage(double& dx, double L) {
    dx -= L * floor(dx / L + 0.5);
}

void periodicWrap(double& x, double L) {
    x -= L * std::nearbyint(x / L);
}