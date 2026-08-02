#include "initializers/xyz_data_loader.h"
#include "units.h"

#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

void XyzDataLoader::loadFromFile(
    const std::string& xyz_filename,
    const std::string& data_unit,
    const std::string& unit_family,
    long init_frame,
    XyzFrameSelectionMode init_frame_mode,
    VecArray& destination,
    const double prefactor
) {
    std::ifstream input_file(xyz_filename);
    if (!input_file.is_open()) {
        throw std::runtime_error(
            std::format("Cannot open the XYZ file named {}.", xyz_filename)
        );
    }

    // XYZ has no frame index, so both selection modes scan complete frames from
    // the beginning of the file.
    for (long frame_index = 0;; ++frame_index) {
        // In step mode, EOF after a complete frame means that the requested
        // step does not exist. EOF inside a frame is reported as malformed.
        if (init_frame_mode == XyzFrameSelectionMode::Step && input_file.peek() == std::char_traits<char>::eof()) {
            throw std::runtime_error(std::format(
                "XYZ file '{}' does not contain a frame with Step {}.",
                xyz_filename,
                init_frame
            ));
        }

        const auto [atom_count, comment] = readFrameHeader(input_file, xyz_filename, frame_index);

        // Index mode compares the ordinal frame. Step mode parses every scanned
        // comment, deliberately enforcing its strict comment-line contract.
        const bool selected = init_frame_mode == XyzFrameSelectionMode::Index
            ? frame_index == init_frame
            : parseStepNumber(xyz_filename, frame_index, comment) == init_frame;

        if (selected) {
            // Only the selected frame is converted to internal units and
            // written into the simulation's data array.
            readFrameBody(
                input_file,
                xyz_filename,
                data_unit,
                unit_family,
                frame_index,
                atom_count,
                destination,
                prefactor
            );
            return;
        }

        // Validate discarded frames as well; raw line skipping could otherwise
        // hide invalid coordinates before the requested frame.
        skipFrameBody(input_file, xyz_filename, frame_index, atom_count);
    }
}

void XyzDataLoader::throwMalformedFrame(
    const std::string& xyz_filename,
    long frame_index,
    const std::string& detail
) {
    throw std::runtime_error(std::format(
        "Malformed XYZ file '{}', frame {}: {}. Expected an atom-count line, "
        "exactly one comment line, and then one data line for each atom.",
        xyz_filename,
        frame_index,
        detail
    ));
}

XyzDataLoader::FrameHeader XyzDataLoader::readFrameHeader(
    std::istream& input,
    const std::string& xyz_filename,
    long frame_index
) {
    std::string atom_count_line;
    if (!std::getline(input, atom_count_line)) {
        throwMalformedFrame(xyz_filename, frame_index, "file ended before the atom-count line");
    }

    if (atom_count_line.empty()) {
        throwMalformedFrame(xyz_filename, frame_index, "found an empty line where the atom-count line was expected");
    }

    std::istringstream atom_count_stream(atom_count_line);
    int atom_count = 0;
    std::string trailing_token;
    if (
        !(atom_count_stream >> atom_count) ||
        atom_count < 0 ||
        atom_count_stream >> trailing_token
        ) {
        throwMalformedFrame(
            xyz_filename,
            frame_index,
            std::format("invalid atom-count line '{}'", atom_count_line)
        );
    }

    std::string comment;
    if (!std::getline(input, comment)) {
        throwMalformedFrame(xyz_filename, frame_index, "file ended before the comment line");
    }

    return FrameHeader{ .atom_count = atom_count, .comment = std::move(comment) };
}

void XyzDataLoader::skipFrameBody(
    std::istream& input,
    const std::string& xyz_filename,
    long frame_index,
    int atom_count
) {
    std::string line;
    for (int atom = 0; atom < atom_count; ++atom) {
        if (!std::getline(input, line)) {
            throwMalformedFrame(
                xyz_filename,
                frame_index,
                std::format("file ended while reading data line {}", atom)
            );
        }

        // Validate skipped records using the same parser as loaded records.
        static_cast<void>(parseDataLine(xyz_filename, frame_index, atom, line));
    }
}

void XyzDataLoader::readFrameBody(
    std::istream& input,
    const std::string& xyz_filename,
    const std::string& data_unit,
    const std::string& unit_family,
    long frame_index,
    int atom_count,
    VecArray& destination,
    double prefactor
) {
    if (destination.len() != atom_count) {
        throw std::runtime_error(std::format(
            "The number of atoms in frame {} of XYZ file '{}' ({}) does not match "
            "the requested number of atoms ({}).",
            frame_index,
            xyz_filename,
            atom_count,
            destination.len()
        ));
    }

    for (int atom = 0; atom < atom_count; ++atom) {
        std::string line;
        if (!std::getline(input, line)) {
            throwMalformedFrame(
                xyz_filename,
                frame_index,
                std::format("file ended while reading data line {}", atom)
            );
        }

        const Vec data = parseDataLine(xyz_filename, frame_index, atom, line);

        for (int axis = 0; axis < NDIM; ++axis) {
            try {
                double converted_value = Units::convertToInternal(
                    unit_family,
                    data_unit,
                    data[axis]
                );

                // For velocities, the prefactor is the mass, which upon multiplication gives the momenta
                if (unit_family == "velocity") {
                    converted_value *= prefactor;
                }

                destination(atom, axis) = converted_value;
            } catch (const std::exception&) {
                throwMalformedFrame(
                    xyz_filename,
                    frame_index,
                    std::format("data line {} cannot be converted to unit '{}': '{}'", atom, data_unit, line)
                );
            }
        }
    }
}

long XyzDataLoader::parseStepNumber(
    const std::string& xyz_filename,
    long frame_index,
    const std::string& comment
) {
    static const std::regex step_pattern(R"(^\s*Step\s+([0-9]+)\s*$)");
    std::smatch match;

    if (!std::regex_match(comment, match, step_pattern)) {
        throwMalformedFrame(
            xyz_filename,
            frame_index,
            std::format(
                "comment line '{}' does not match the required 'Step <non-negative integer>' format for step selection",
                comment
            )
        );
    }

    try {
        return std::stol(match[1].str());
    } catch (const std::exception&) {
        throwMalformedFrame(
            xyz_filename,
            frame_index,
            std::format("step number '{}' is outside the supported range", match[1].str())
        );
    }
}

Vec XyzDataLoader::parseDataLine(
    const std::string& xyz_filename,
    long frame_index,
    int atom_index,
    const std::string& line
) {
    std::istringstream line_stream(line);
    std::vector<std::string> tokens;
    std::string token;
    while (line_stream >> token) {
        tokens.push_back(token);
    }

    // XYZ file line should have exactly one atom label, followed by 3 numeric values (even when NDIM is 2 or 1, the file format still expects 3 values for compatibility)
    if (tokens.size() < 4) {
        throwMalformedFrame(
            xyz_filename,
            frame_index,
            std::format(
                "data line {} is '{}'; expected an atom identifier followed by {} values",
                atom_index,
                line,
                NDIM
            )
        );
    }

    Vec data{};
    for (int axis = 0; axis < NDIM; ++axis) {
        const std::string& data_token = tokens[1 + axis];  // Always start from position 1 (after atom label)

        try {
            std::size_t parsed_length = 0;
            data[axis] = std::stod(data_token, &parsed_length);
            if (parsed_length != data_token.size()) {
                throw std::invalid_argument("trailing characters");
            }
        } catch (const std::exception&) {
            throwMalformedFrame(
                xyz_filename,
                frame_index,
                std::format("data line {} contains an invalid value: '{}'", atom_index, line)
            );
        }
    }

    return data;
}