#pragma once

#include "core/simulation_config.h"

#include <istream>
#include <string>

/**
 * @brief Utility for loading data from XYZ files with proper validation and unit conversion.
 *
 * This utility generalizes the XYZ file parsing logic to support both position and velocity data.
 * The primary difference is in the unit conversion applied (length vs velocity).
 */
class XyzDataLoader {
public:
    /**
     * The two non-data records that start every strict XYZ frame.
     *
     * The comment is retained because step-based selection inspects it before
     * deciding whether to load or skip the frame body.
     */
    struct FrameHeader {
        int atom_count;
        std::string comment;
    };

    /**
     * Loads data from an XYZ file with specified unit conversion.
     *
     * @param xyz_filename Path to the XYZ file
     * @param data_unit Unit declared in the configuration for XYZ data values
     * @param unit_family The unit family for conversion ("length" or "velocity")
     * @param init_frame Frame index or step number to load
     * @param init_frame_mode Determines how init_frame is interpreted (Index or Step)
     * @param destination Target array to store loaded data
     * @param prefactor Constant pre-factor applied to all loaded values (default is 1.0)
     */
    static void loadFromFile(
        const std::string& xyz_filename,
        const std::string& data_unit,
        const std::string& unit_family,
        long init_frame,
        XyzFrameSelectionMode init_frame_mode,
        VecArray& destination,
        const double prefactor = 1.0
    );

private:
    [[noreturn]] static void throwMalformedFrame(
        const std::string& xyz_filename,
        long frame_index,
        const std::string& detail
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
     * Consumes and validates data records without storing them.
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
        const std::string& data_unit,
        const std::string& unit_family,
        long frame_index,
        int atom_count,
        VecArray& destination,
        double prefactor
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
     * Validates one atom record and returns its final NDIM data fields.
     *
     * This preserves the existing convention: leading atom labels and metadata
     * are accepted, while the final NDIM fields hold the data values.
     */
    static Vec parseDataLine(
        const std::string& xyz_filename,
        long frame_index,
        int atom_index,
        const std::string& line
    );
};