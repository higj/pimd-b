#include "initializers/xyz_position_initializer.h"

#include <fstream>

XyzPositionInitializer::XyzPositionInitializer(
    const std::string& filename, int first_idx, const std::shared_ptr<dVec>& coord, double box_size)
    : PositionInitializer(coord, box_size), m_filename(filename), m_first_idx(first_idx) {
}

void XyzPositionInitializer::initialize() {
    //const int arg = m_state->currentBead() + m_first_idx;
    //loadFromFile(std::vformat(m_filename, std::make_format_args(arg)), m_coord);
    loadFromFile(std::vformat(m_filename, std::make_format_args(m_first_idx)), *m_coord);
}

void XyzPositionInitializer::loadFromFile(const std::string& pos_filename, dVec& destination) const {
    std::ifstream input_file(m_filename);

    if (!input_file.is_open())
        throw std::runtime_error(std::format("Cannot open the xyz file named {}.", m_filename));

    int numAtoms;
    input_file >> numAtoms;

    if (destination.len() != numAtoms)
        throw std::runtime_error(
            std::format("The number of atoms in the xyz file ({}) does not match the requested number of atoms.",
                m_filename)
        );

    // The command "inputFile >> numAtoms" read the number, but the newline character that follows that number was not consumed. 
    // Later, we use std::getline(inputFile, comment) to read the comment line, so to prevent it from
    // stopping immediately (due to encountering the newline character from the previous line), we must consume
    // the newline character by using ignore (the first argument ensures that we read as many characters as needed until we
    // find the delimiter '\n' for the first time).
    input_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Currently we ignore the comment line.
    // TODO: Extract the units or the cell parameters from the comment line
    // TODO: If the provided xyz file does not contain a comment line, the first coordinate
    // will be treated as a comment and therefore will be lost. Fix this.
    std::string comment;

    std::getline(input_file, comment, '\n');

    for (int i = 0; i < numAtoms; ++i) {
        std::string symbol;

        // Currently we ignore the symbol
        input_file >> symbol;

        // TODO: Add a check that the number of dimensions in the .xyz file matches NDIM
        for (int j = 0; j < NDIM; ++j) {
            input_file >> destination(i, j);
            destination(i, j) = Units::convertToInternal("length", "angstrom", destination(i, j));
        }
    }

    input_file.close();
}