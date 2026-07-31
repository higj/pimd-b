#include "initializers/xyz_position_initializer.h"

#include "units.h"

#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

// Construction is intentionally side-effect free; file I/O is deferred until
// initialize(), after the simulation has created its state and contexts.
XyzPositionInitializer::XyzPositionInitializer(
    const std::string& filename,
    int first_idx,
    const std::string& init_pos_unit,
    long init_pos_frame,
    XyzFrameSelectionMode init_pos_frame_mode,
    const std::shared_ptr<VecArray>& coord,
    const BoxContext& box_ctx
)
    : PositionInitializer(coord, box_ctx),
      m_filename(filename),
      m_first_idx(first_idx),
      m_init_pos_unit(init_pos_unit),
      m_init_pos_frame(init_pos_frame),
      m_init_pos_frame_mode(init_pos_frame_mode) {
}

void XyzPositionInitializer::initialize() {
    // A literal filename is valid: vformat returns it unchanged when there are
    // no replacement fields. A formatted name supports per-bead XYZ files.
    const std::string xyz_filename = std::vformat(
        m_filename,
        std::make_format_args(m_first_idx)
    );

    loadFromFile(
        xyz_filename,
        m_init_pos_unit,
        m_init_pos_frame,
        m_init_pos_frame_mode,
        *m_coord
    );
}

void XyzPositionInitializer::throwMalformedFrame(
    const std::string& xyz_filename,
    long frame_index,
    const std::string& detail
) {
    // Use one message shape for all syntax failures, so callers and MPI error
    // handling always retain the physical file and zero-based frame location.
    throw std::runtime_error(std::format(
        "Malformed XYZ file '{}', frame {}: {}. Expected an atom-count line, "
        "exactly one comment line, and then one coordinate line for each atom.",
        xyz_filename,
        frame_index,
        detail
    ));
}

void XyzPositionInitializer::loadFromFile(
    const std::string& xyz_filename,
    const std::string& init_pos_unit,
    long init_pos_frame,
    XyzFrameSelectionMode init_pos_frame_mode,
    VecArray& destination
) {
    // Params validates that init_pos_frame is non-negative before construction.
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
        if (init_pos_frame_mode == XyzFrameSelectionMode::Step && input_file.peek() == std::char_traits<char>::eof()) {
            throw std::runtime_error(std::format(
                "XYZ file '{}' does not contain a frame with Step {}.",
                xyz_filename,
                init_pos_frame
            ));
        }

        const auto [atom_count, comment] = readFrameHeader(input_file, xyz_filename, frame_index);

        // Index mode compares the ordinal frame. Step mode parses every scanned
        // comment, deliberately enforcing its strict comment-line contract.
        const bool selected = init_pos_frame_mode == XyzFrameSelectionMode::Index
            ? frame_index == init_pos_frame
            : parseStepNumber(xyz_filename, frame_index, comment) == init_pos_frame;

        if (selected) {
            // Only the selected frame is converted to internal units and
            // written into the simulation's coordinate array.
            readFrameBody(
                input_file,
                xyz_filename,
                init_pos_unit,
                frame_index,
                atom_count,
                destination
            );
            return;
        }

        // Validate discarded frames as well; raw line skipping could otherwise
        // hide invalid coordinates before the requested frame.
        skipFrameBody(input_file, xyz_filename, frame_index, atom_count);
    }
}

XyzPositionInitializer::FrameHeader XyzPositionInitializer::readFrameHeader(
    std::istream& input,
    const std::string& xyz_filename,
    long frame_index
) {
    // getline is used throughout to make unexpected blank lines observable.
    // Formatted extraction would silently skip whitespace and mask the error.
    std::string atom_count_line;
    if (!std::getline(input, atom_count_line)) {
        throwMalformedFrame(xyz_filename, frame_index, "file ended before the atom-count line");
    }

    if (atom_count_line.empty()) {
        throwMalformedFrame(xyz_filename, frame_index, "found an empty line where the atom-count line was expected");
    }

    // A valid atom-count record contains exactly one non-negative integer.
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

    // Standard XYZ includes exactly one comment line, which may be empty.
    std::string comment;
    if (!std::getline(input, comment)) {
        throwMalformedFrame(xyz_filename, frame_index, "file ended before the comment line");
    }

    return FrameHeader{ .atom_count = atom_count, .comment = std::move(comment) };
}

void XyzPositionInitializer::skipFrameBody(
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
                std::format("file ended while reading coordinate line {}", atom)
            );
        }

        // Validate skipped records using the same parser as loaded records.
        static_cast<void>(parseCoordinateLine(xyz_filename, frame_index, atom, line));
    }
}

void XyzPositionInitializer::readFrameBody(
    std::istream& input,
    const std::string& xyz_filename,
    const std::string& init_pos_unit,
    long frame_index,
    int atom_count,
    VecArray& destination
) {
    // The simulation uses a fixed particle count, so a differently sized frame
    // cannot be stored safely in destination.
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
                std::format("file ended while reading coordinate line {}", atom)
            );
        }

        // Parse before converting so a syntax failure is distinct from a unit
        // conversion failure in the diagnostic message.
        const Vec coordinates = parseCoordinateLine(xyz_filename, frame_index, atom, line);

        for (int axis = 0; axis < NDIM; ++axis) {
            try {
                destination(atom, axis) = Units::convertToInternal(
                    "length",
                    init_pos_unit,
                    coordinates[axis]
                );
            } catch (const std::exception&) {
                throwMalformedFrame(
                    xyz_filename,
                    frame_index,
                    std::format("coordinate line {} cannot be converted to unit '{}': '{}'", atom, init_pos_unit, line)
                );
            }
        }
    }
}

long XyzPositionInitializer::parseStepNumber(
    const std::string& xyz_filename,
    long frame_index,
    const std::string& comment
) {
    // The anchors reject comments such as "Before Step 10" or "Step 10 done".
    // Step selection is therefore deterministic rather than a substring match.
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
        // stol also detects a syntactically valid value that cannot fit in long.
        return std::stol(match[1].str());
    } catch (const std::exception&) {
        throwMalformedFrame(
            xyz_filename,
            frame_index,
            std::format("step number '{}' is outside the supported range", match[1].str())
        );
    }
}

Vec XyzPositionInitializer::parseCoordinateLine(
    const std::string& xyz_filename,
    long frame_index,
    int atom_index,
    const std::string& line
) {
    // XYZ variants often include labels or metadata before coordinates. The
    // existing input convention uses the final NDIM tokens as coordinates.
    std::istringstream line_stream(line);
    std::vector<std::string> tokens;
    std::string token;
    while (line_stream >> token) {
        tokens.push_back(token);
    }

    if (tokens.size() < NDIM + 1) {
        throwMalformedFrame(
            xyz_filename,
            frame_index,
            std::format(
                "coordinate line {} is '{}'; expected an atom identifier followed by {} coordinate values",
                atom_index,
                line,
                NDIM
            )
        );
    }

    Vec coordinates{};
    for (int axis = 0; axis < NDIM; ++axis) {
        const std::string& coordinate_token = tokens[tokens.size() - NDIM + axis];

        try {
            // stod accepts numeric prefixes such as "1.0abc". Requiring full
            // consumption prevents those malformed fields from being accepted.
            std::size_t parsed_length = 0;
            coordinates[axis] = std::stod(coordinate_token, &parsed_length);
            if (parsed_length != coordinate_token.size()) {
                throw std::invalid_argument("trailing characters");
            }
        } catch (const std::exception&) {
            throwMalformedFrame(
                xyz_filename,
                frame_index,
                std::format("coordinate line {} contains an invalid value: '{}'", atom_index, line)
            );
        }
    }

    return coordinates;
}
