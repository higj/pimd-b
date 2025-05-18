#include "initializers/manual_momentum_initializer.h"

#include <fstream>

ManualMomentumInitializer::ManualMomentumInitializer(
    const std::string& filename,
    int first_idx,
    const std::shared_ptr<SystemState>& state,
    double mass) : MomentumInitializer(state, mass), m_filename(filename), m_first_idx(first_idx)
{
}

void ManualMomentumInitializer::initialize()
{
    //const std::string vel_filename_format = std::get<std::string>(sim_params.at("init_vel_manual_filename_format"));
    //const int arg = m_state->currentBead() + std::get<int>(sim_params.at("init_vel_first_index"));
    const int arg = m_state->currentBead() + m_first_idx;

    loadFromFile(std::vformat(m_filename, std::make_format_args(arg)), m_state->momenta);
}

void ManualMomentumInitializer::loadFromFile(const std::string& vel_filename, dVec& destination) const
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

    for (int i = 0; i < num_atoms; ++i) {
        std::string symbol;

        // We consume twice, because typical LAMMPS velocity dumps contain lines such as:
        // ATOM_NUM ID v_x v_y v_z
        input_file >> symbol;
        input_file >> symbol;

        // TODO: Add a check that the number of dimensions in the .xyz file matches NDIM
        for (int j = 0; j < NDIM; ++j) {
            input_file >> destination(i, j);
            // For LAMMPS velocity files
            destination(i, j) = m_mass * Units::convertToInternal("velocity", "angstrom/ps", destination(i, j));

            // For i-Pi velocity files
            //destination(i, j) = mass * destination(i, j);
        }
    }

    input_file.close();
}
