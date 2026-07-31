#pragma once

#include "core/simulation_config.h"
#include "initializers/position_initializer.h"

#include <istream>
#include <memory>
#include <string>

class XyzPositionInitializer final : public PositionInitializer {
public:
    /**
     * @brief Creates an initializer that reads one frame from an XYZ trajectory.
     *
     * The filename may have one format replacement field, which is replaced
     * with first_idx so individual MPI beads can use separate input files.
     * init_pos_frame is either a zero-based frame index or a value from a
     * "Step <number>" comment, according to init_pos_frame_mode.
     */
    XyzPositionInitializer(
        const std::string& filename,
        int first_idx,
        const std::string& init_pos_unit,
        long init_pos_frame,
        XyzFrameSelectionMode init_pos_frame_mode,
        const std::shared_ptr<VecArray>& coord,
        const BoxContext& box_ctx
    );

    /**
     * Initializes positions from user-provided files.
     */
    void initialize() override;

private:
    std::string m_filename;
    // Index used to format m_filename for the current MPI bead.
    int m_first_idx;
    // Unit declared in the configuration for XYZ coordinate values.
    std::string m_init_pos_unit;
    // User-requested frame index or step number.
    long m_init_pos_frame;
    // Determines how m_init_pos_frame is interpreted.
    XyzFrameSelectionMode m_init_pos_frame_mode;

    /**
     * The two non-coordinate records that start every strict XYZ frame.
     *
     * The comment is retained because step-based selection inspects it before
     * deciding whether to load or skip the frame body.
     */
    struct FrameHeader {
        int atom_count;
        std::string comment;
    };

    [[noreturn]] static void throwMalformedFrame(
        const std::string& xyz_filename,
        long frame_index,
        const std::string& detail
    );

    /**
     * Loads the requested frame from an XYZ file.
     */
    static void loadFromFile(
        const std::string& xyz_filename,
        const std::string& init_pos_unit,
        long init_pos_frame,
        XyzFrameSelectionMode init_pos_frame_mode,
        VecArray& destination
    );

    /**
     * Reads and validates the atom-count and single comment line of a frame.
     */
    static FrameHeader readFrameHeader(
        std::istream& input,
        const std::string& xyz_filename,
        long frame_index
    );

    /**
     * Consumes and validates coordinate records without storing them.
     *
     * Skipped frames are deliberately validated so malformed earlier frames
     * cannot silently change the interpretation of a requested later frame.
     */
    static void skipFrameBody(
        std::istream& input,
        const std::string& xyz_filename,
        long frame_index,
        int atom_count
    );

    /**
     * Reads, validates, converts, and stores all records in the selected frame.
     */
    static void readFrameBody(
        std::istream& input,
        const std::string& xyz_filename,
        const std::string& init_pos_unit,
        long frame_index,
        int atom_count,
        VecArray& destination
    );

    /**
     * Extracts a numeric value from a strict "Step <non-negative integer>" comment.
     */
    static long parseStepNumber(
        const std::string& xyz_filename,
        long frame_index,
        const std::string& comment
    );

    /**
     * Validates one atom record and returns its final NDIM coordinate fields.
     *
     * This preserves the existing convention: leading atom labels and metadata
     * are accepted, while the final NDIM fields hold the coordinates.
     */
    static Vec parseCoordinateLine(
        const std::string& xyz_filename,
        long frame_index,
        int atom_index,
        const std::string& line
    );
};
