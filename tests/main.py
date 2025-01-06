from math import exp
import re
import subprocess
import os
import shutil
from pathlib import Path
import sys
import numpy as np
import configparser
import argparse
import MDAnalysis as mda


# Columns of observables to compare
#columns = ["step", "kinetic", "potential", "ext_pot", "int_pot", "virial", "temperature", "cl_kinetic", "cl_spring"]
out_filename = "simulation.out"


def read_data(file_path):
    # Read the data from the output file
    data = np.genfromtxt(file_path, names=True, skip_header=0)
    column_names = data.dtype.names

    # Extract columns and values
    values = [data[i] for i in column_names]

    # Create a dictionary to store the data
    data_dict = {column: values[i] for i, column in enumerate(column_names)}

    return data_dict

def get_number_of_beads(input_file):
    config = configparser.ConfigParser()
    
    try:
        config.read(input_file)
        
        # Extract the value associated with the "nbeads" key
        nbeads_value = int(float(config.get('simulation', 'nbeads')))
        
        return nbeads_value
    except configparser.MissingSectionHeaderError:
        raise RuntimeError(f"Test '{input_file}' does not contain the [simulation] section!")
    except configparser.NoOptionError:
        raise RuntimeError(f"Test '{input_file}' does not contain the 'nbeads' key in the [simulation] section!")


def get_number_of_atoms(input_file):
    config = configparser.ConfigParser()
    
    try:
        config.read(input_file)
        
        # Extract the value associated with the "nbeads" key
        nbeads_value = int(float(config.get('system', 'natoms')))
        
        return nbeads_value
    except configparser.MissingSectionHeaderError:
        raise RuntimeError(f"Test '{input_file}' does not contain the [system] section!")
    except configparser.NoOptionError:
        raise RuntimeError(f"Test '{input_file}' does not contain the 'natoms' key in the [system] section!")


def run_simulation(executable_dir, input_file):
    nbeads = get_number_of_beads(input_file)
    
    mpirun_bin_dir = os.environ['OPENMPI_BIN']
    mpi = Path(mpirun_bin_dir) / 'mpirun'
    mpirun_command = [mpi, '--oversubscribe', '-np', str(nbeads), (executable_dir / 'pimdb').as_posix(), '-partition', f"{nbeads}x1", '-in', input_file.as_posix()]
    
    # Windows version
    #mpirun_command = ['mpiexec', '-n', str(nbeads), (executable_dir / 'pimdb.exe').as_posix(), '-in', input_file.as_posix()]
    
    return subprocess.run(mpirun_command, capture_output=True, text=True, check=True)


def compare_arrays(array1, array2):
    # Check if the arrays are equal within a tolerance
    if np.allclose(array1, array2, rtol=1e-5):
        return True, None  # Arrays are equal

    # Find the index where the first difference occurs
    index = np.where(~np.isclose(array1, array2))[0][0]

    # Return the difference
    return False, index


def compare_output(actual_output, expected_output):
    data_actual = read_data(actual_output)
    data_expected = read_data(expected_output)
    
    # Check if the number of steps match
    if len(data_actual['step']) != len(data_expected['step']):
        raise AssertionError("Test failed: Number of steps do not match.")
    
    columns = list(data_actual.keys())
    
    print(f"Comparing columns {columns}")
    
    # Check if the values of observables match
    for column in columns:
        are_equal, index = compare_arrays(data_actual[column], data_expected[column])
        if not are_equal:
            print(f"Test 'failed': '{column}' does not match at step {index}.")
           # raise AssertionError(f"Test failed: '{column}' does not match at step {index}.")
        
    return True


def get_coordinates(xyz_path, natoms, dim=3):
    xyz_class = mda.coordinates.XYZ.XYZReader(xyz_path)
    steps = xyz_class.n_frames

    coords = np.zeros((steps, natoms, dim))

    for timestep in xyz_class:
        coords[timestep.frame, :, :] = timestep.positions

    return coords, steps


# Coordinates are assumed to be in Angstroms
# Coordinate files are assumed to have names in the format "position_0.xyz", ..., "position_P-1.xyz"
def compare_xyz(actual_xyz_file, expected_xyz_file, natoms):
    coords, steps = get_coordinates(xyz_path=actual_xyz_file, natoms=natoms, dim=3)
    coords_test, steps_test = get_coordinates(xyz_path=expected_xyz_file, natoms=natoms, dim=3)

    are_equal, index = compare_arrays(coords, coords_test)
    if not are_equal:
        print(coords[index:],coords_test[index:])
        raise AssertionError(f"Test failed: Coordinates in {actual_xyz_file.name} do not match at step {index}.")
    
    return True


def is_float(element: any) -> bool:
    if element is None:
        return False
    try:
        float(element)
        return True
    except ValueError:
        return False


def extract_numeric_data(file_path):
    """
    Extract numeric data from the file, focusing on the last 3 values in each row.
    """
    numeric_data = []

    with open(file_path, 'r') as file:
        for line_idx, line in enumerate(file):
            values = line.strip().split()[-3:]
            vals_are_floats = [is_float(value) for value in values]
            
            # Check if the values are valid numbers
            # Also check that the velocity vector has at least three components
            if all(vals_are_floats) and len(vals_are_floats) > 2:
                numeric_data.append([float(value) for value in values])

    return np.array(numeric_data)


def test_coordinates(output_folder, test_folder, in_file):
    test_xyz_files = list(test_folder.glob("*.xyz"))
            
    # If the test contains '.xyz' files, compare them.
    # Otherwise, skip this test.
    if test_xyz_files:
        print("Comparing trajectories...")
        xyz_file_names = [file.name for file in test_xyz_files]  # Expected xyz file names

        # Check if the simulation output directory has the same number of xyz files and same filenames
        out_xyz_files = list(output_folder.glob("*.xyz"))
                
        if len(out_xyz_files) != len(test_xyz_files):
            raise AssertionError(f"Test failed: number of xyz files found in {output_folder} is incorrect.")
                
        for xyz_file in out_xyz_files:
            if xyz_file.name not in xyz_file_names:
                raise AssertionError(f"Test failed: the generated xyz files have incorrect names.")
                
        natoms = get_number_of_atoms(tests_dir / in_file)
                
        for xyz_file_name in xyz_file_names:
            out_xyz_file = output_folder / xyz_file_name
            test_xyz_file = test_folder / xyz_file_name
            compare_xyz(actual_xyz_file=out_xyz_file, expected_xyz_file=test_xyz_file, natoms=natoms)
                    
        print("Test passed: Trajectories match.")


def test_velocities(output_folder, test_folder):
    # Function to filter files with names in the format "velocity_X.dat"
    def is_velocity_file(file):
        return file.name.startswith("velocity_") and file.name.endswith(".dat") and file.stem.split("_")[1].isdigit()

    # Check if the test case has at least one file with the format "velocity_X.dat"
    velocity_files = list(filter(is_velocity_file, test_folder.glob("velocity_*.dat")))
    
    # If not, quit the test
    if not velocity_files:
        return True

    print("Comparing velocities...")
    
    # Get the names of all 'velocity_X.dat' files in test
    velocity_file_names = [file.name for file in velocity_files]

    # Check if the output directory has the same number of files and same filenames
    out_velocity_files = list(filter(is_velocity_file, output_folder.glob("velocity_*.dat")))
    if len(out_velocity_files) != len(velocity_files):
        raise AssertionError(f"Test failed: Different number of 'velocity_X.dat' files found in {output_folder}")
    
    for velocity_file in out_velocity_files:
        if velocity_file.name not in velocity_file_names:
            raise AssertionError(f"Test failed: Different filenames found in {output_folder}")

    # Iterate over similar files in both directories and compare their contents
    for velocity_file_name in velocity_file_names:
        test_vel_file = test_folder / velocity_file_name
        out_vel_file = output_folder / velocity_file_name
        
        vels_out = extract_numeric_data(out_vel_file)
        vels_test = extract_numeric_data(test_vel_file)

        are_equal, index = compare_arrays(vels_out, vels_test)
        if not are_equal:
            raise AssertionError(f"Test failed: Velocities do not match at step {index}.")
    
    print("Test passed: Velocities match.")

    return True

def test_forces(output_folder, test_folder):
    # Function to filter files with names in the format "force_X.dat"
    def is_force_file(file):
        return file.name.startswith("force_") and file.name.endswith(".dat") and file.stem.split("_")[1].isdigit()

    # Check if the test case has at least one file with the format "force_X.dat"
    force_files = list(filter(is_force_file, test_folder.glob("force_*.dat")))
    
    # If not, quit the test
    if not force_files:
        return True

    print("Comparing forces...")
    
    # Get the names of all 'force_X.dat' files in test
    force_file_names = [file.name for file in force_files]

    # Check if the output directory has the same number of files and same filenames
    out_force_files = list(filter(is_force_file, output_folder.glob("force_*.dat")))
    if len(out_force_files) != len(force_files):
        raise AssertionError(f"Test failed: Different number of 'force_X.dat' files found in {output_folder}")
    
    for force_file in out_force_files:
        if force_file.name not in force_file_names:
            raise AssertionError(f"Test failed: Different filenames found in {output_folder}")

    # Iterate over similar files in both directories and compare their contents
    for force_file_name in force_file_names:
        test_force_file = test_folder / force_file_name
        out_force_file = output_folder / force_file_name
        
        forces_out = extract_numeric_data(out_force_file)
        forces_test = extract_numeric_data(test_force_file)

        are_equal, index = compare_arrays(forces_out, forces_test)
        if not are_equal:
            raise AssertionError(f"Test failed: Forces do not match at step {index}.")
    
    print("Test passed: Forces match.")

    return True


def run_tests(executable_dir, tests_dir, is_old_bosonic):
    # Navigate to the directory containing test cases
    os.chdir(tests_dir)
    
    # The output folder will be generated in the current working directory
    out_folder = Path.cwd() / 'output'
    out_path = out_folder / out_filename
    
    # If old_bosonic is True, iterate over any test case whose folder does not have the prefix "bosonic_quadratic".
    # If old_bosonic is False, iterate over any test case whose folder does not have the prefix "bosonic_factorial".
    
    excluded_bosonic_prefix = "bosonic_quadratic" if is_old_bosonic else "bosonic_factorial"
    test_cases = [test_case for test_case in tests_dir.iterdir() if test_case.is_dir() and not test_case.name.startswith(excluded_bosonic_prefix)]
    
    # Iterate over each test case folder
    for test_case in test_cases:
        #print(test_case)
        print("-------------------")
        #print("Currently here:", os.getcwd())
        if os.path.isdir(test_case):
            print(f"Running test case: {test_case.name}")
            input_file = test_case / f"{test_case.name}.ini"
            # print("Input file of the test:", input_file)
            expected_out_path = test_case / out_filename
            
            # Run the simulation and get the path to the output file
            cmd_out = run_simulation(executable_dir, tests_dir / input_file)
            #print(cmd_out)
            
            print("Comparing output files...")
            # 1st test: Compare the output with the expected output
            if compare_output(out_path, expected_out_path):
                print("Test passed: Output matches expected output.")
            else:
                raise AssertionError("Test failed: Output does not match expected output.")
            
            # 2nd test: Compare the coordinates with the expected coordinates
            test_coordinates(out_folder, test_case, input_file)
            
            # 3rd test: Compare the velocities with the expected velocities
            test_velocities(out_folder, test_case)
            
            # 4th test: Compare the forces with the expected forces
            test_forces(out_folder, test_case)
            
            # Clean up the generated output file
            print("Deleting:", out_folder)
            shutil.rmtree(out_folder)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PIMD-B++ tests')
    parser.add_argument('executable_dir', type=Path, help='Directory containing the executable')
    parser.add_argument('tests_dir', type=Path, help='Directory containing the test cases')
    parser.add_argument('--old-bosonic', default=False, action=argparse.BooleanOptionalAction, 
                        help='Mark if the executable was compiled with the old bosonic algorithm')

    args = parser.parse_args()

    executable_dir = Path(args.executable_dir)
    tests_dir = Path(args.tests_dir)
    is_old_bosonic = args.old_bosonic
    
    # Run tests
    run_tests(executable_dir, tests_dir, is_old_bosonic)
