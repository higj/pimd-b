#include <iostream>
#include <fstream>

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

void loadTrajectories(const std::string& xyz_filename, dVec& destination) {
    std::ifstream inputFile(xyz_filename);

    if (!inputFile.is_open())
        throw std::runtime_error(std::format("Cannot open the xyz file named {}.", xyz_filename));

    int numAtoms;
    inputFile >> numAtoms;

    if (destination.len() != numAtoms)
        throw std::runtime_error(
            std::format("The number of atoms in the xyz file ({}) does not match the requested number of atoms.",
                        xyz_filename)
        );

    // The command "inputFile >> numAtoms" read the number, but the newline character that follows that number was not consumed. 
    // Later, we use std::getline(inputFile, comment) to read the comment line, so to prevent it from
    // stopping immediately (due to encountering the newline character from the previous line), we must consume
    // the newline character by using ignore (the first argument ensures that we read as many characters as needed until we
    // find the delimiter '\n' for the first time).
    inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Currently we ignore the comment line.
    // TODO: Extract the units or the cell parameters from the comment line
    // TODO: If the provided xyz file does not contain a comment line, the first coordinate
    // will be treated as a comment and therefore will be lost. Fix this.
    std::string comment;

    std::getline(inputFile, comment, '\n');

    for (int i = 0; i < numAtoms; ++i) {
        std::string symbol;

        // Currently we ignore the symbol
        inputFile >> symbol;

        // TODO: Add a check that the number of dimensions in the .xyz file matches NDIM
        for (int j = 0; j < NDIM; ++j) {
            inputFile >> destination(i, j);
            destination(i, j) = Units::convertToInternal("length", "angstrom", destination(i, j));
        }
    }

    inputFile.close();
}

void loadMomenta(const std::string& vel_filename, double mass, dVec& destination) {
    std::ifstream inputFile(vel_filename);

    if (!inputFile.is_open())
        throw std::runtime_error(std::format("Cannot open the velocity file named {}.", vel_filename));

    int numAtoms;
    inputFile >> numAtoms;

    if (destination.len() != numAtoms)
        throw std::runtime_error(
            std::format("The number of atoms in the velocity file ({}) does not match the requested number of atoms.",
                        vel_filename)
        );

    inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string comment;

    std::getline(inputFile, comment, '\n');

    for (int i = 0; i < numAtoms; ++i) {
        std::string symbol;

        // We consume twice, because typical LAMMPS velocity dumps contain lines such as:
        // ATOM_NUM ID v_x v_y v_z
        inputFile >> symbol;
        inputFile >> symbol;

        // TODO: Add a check that the number of dimensions in the .xyz file matches NDIM
        for (int j = 0; j < NDIM; ++j) {
            inputFile >> destination(i, j);
            // For LAMMPS velocity files
            destination(i, j) = mass * Units::convertToInternal("velocity", "angstrom/ps", destination(i, j));

            // For i-Pi velocity files
            //destination(i, j) = mass * destination(i, j);
        }
    }

    inputFile.close();
}
