#include "initializers/manual_momentum_initializer.h"

#include <fstream>
#include <sstream>

ManualMomentumInitializer::ManualMomentumInitializer(
    const std::string& filename,
    int first_idx,
    const std::string& init_vel_unit,
    const std::shared_ptr<SystemState>& state,
    double mass) : MomentumInitializer(state, mass), m_filename(filename), m_first_idx(first_idx), m_init_vel_unit(init_vel_unit)
{
}

void ManualMomentumInitializer::initialize()
{
    //const std::string vel_filename_format = std::get<std::string>(sim_params.at("init_vel_manual_filename_format"));
    //const int arg = m_state->currentBead() + std::get<int>(sim_params.at("init_vel_first_index"));
    const int arg = m_state->currentBead() + m_first_idx;

    loadFromFile(std::vformat(m_filename, std::make_format_args(arg)), m_init_vel_unit, m_state->momenta);
}

void ManualMomentumInitializer::loadFromFile(const std::string& vel_filename, const std::string& init_vel_unit, VecArray& destination) const
{
    std::ifstream input_file(vel_filename);

    if (!input_file.is_open())
        throw std::runtime_error(std::format("Cannot open the velocity file named {}.", vel_filename));

    int num_atoms;
    input_file >> num_atoms;

    if (destination.len() != num_atoms)
        throw std::runtime_error(
            std::format("The number of atoms in the velocity file ({}) does not match the requested number of atoms.",
                vel_filename)
        );

    input_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string comment;

    std::getline(input_file, comment, '\n');

    /*
    for (int i = 0; i < num_atoms; ++i) {
        std::string symbol;

        // We consume twice, because typical LAMMPS velocity dumps contain lines such as:
        // ATOM_NUM ID v_x v_y v_z
        input_file >> symbol;
        input_file >> symbol;

        // TODO: Add a check that the number of dimensions in the .xyz file matches NDIM
        for (int j = 0; j < NDIM; ++j) {
            input_file >> destination(i, j);
            // Note: LAMMPS velocity files use "angstrom/ps", i-PI velocity files use "atomic_units"
            destination(i, j) = m_mass * Units::convertToInternal("velocity", init_vel_unit, destination(i, j));
        }

        // Skip the rest of the line (remaining velocities if NDIM < 3)
        input_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    */
    for (int i = 0; i < num_atoms; ++i) {
        std::string line;
        if (!std::getline(input_file, line)) {
            throw std::runtime_error("Unexpected end of file while reading velocities");
        }

        // Split the line into tokens
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        while (iss >> token) {
            tokens.push_back(token);
        }

        // Validate we have at least NDIM + 1 elements (at least one ID + velocities)
        if (tokens.size() < NDIM + 1) {
            throw std::runtime_error("Line " + std::to_string(i + 1) +
                " has insufficient fields: expected at least " +
                std::to_string(NDIM + 1) + ", got " +
                std::to_string(tokens.size()));
        }

        // Extract the last NDIM values as velocities
        for (int j = 0; j < NDIM; ++j) {
            try {
                destination(i, j) = std::stod(tokens[tokens.size() - NDIM + j]);
                // Note: LAMMPS velocity files use "angstrom/ps", i-PI velocity files use "atomic_units"
                destination(i, j) = m_mass * Units::convertToInternal("velocity", init_vel_unit, destination(i, j));
            } catch ([[maybe_unused]] const std::invalid_argument& e) {
                throw std::runtime_error("Invalid velocity value at line " +
                    std::to_string(i + 1) + ", component " +
                    std::to_string(j));
            }
        }
    }

    input_file.close();
}
