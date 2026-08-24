#include <iostream>

#include "common.h"

void printMsg(const std::string& msg, int this_rank, int out_rank) {
    if (this_rank == out_rank) {
        std::cout << msg << '\n';
    }
}

void printWarning(const std::string& msg, int this_rank, int out_rank) {
    //printMsg("[!] " + msg, this_rank, out_rank);
    if (this_rank == out_rank) {
        std::cerr << "[!] Warning: " << msg << '\n';
    }
}

/*
void printInfo(const std::string& info, bool& info_flag, int this_rank, int out_rank) {
    printMsg(info, this_rank, out_rank);
    info_flag = true;
}
*/

void printStatus(const std::string& status, int this_rank, int out_rank) {
    printMsg("[*] " + status, this_rank, out_rank);
}

void printError(const std::string& msg, int this_rank, const std::string& err_type, int out_rank) {
    printMsg("[X] " + err_type + ": " + msg, this_rank, out_rank);
}
