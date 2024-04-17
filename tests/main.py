import subprocess
import os
import shutil
from pathlib import Path
import sys
import numpy as np
import configparser
import argparse


# Columns of observables to compare
columns = ["step", "kinetic", "potential", "ext_pot", "int_pot", "virial"]
out_filename = "simulation.out"


def read_data(file_path):
    # Read the data from the output file
    data = np.genfromtxt(file_path, skip_header=1)

    # Extract columns and values
    values = [data[:, i] for i in range(len(columns))]

    # Create a dictionary to store the data
    data_dict = {column: values[i] for i, column in enumerate(columns)}

    return data_dict


def get_number_of_beads(input_file):
    config = configparser.ConfigParser()
    
    try:
        config.read(input_file)
        
        # Extract the value associated with the "nbeads" key
        nbeads_value = int(config.get('simulation', 'nbeads'))
        
        return nbeads_value
    except configparser.MissingSectionHeaderError:
        raise RuntimeError(f"Test '{input_file}' does not contain the [simulation] section!")
    except configparser.NoOptionError:
        raise RuntimeError(f"Test '{input_file}' does not contain the 'nbeads' key in the [simulation] section!")


def run_simulation(executable_dir, input_file):
    nbeads = get_number_of_beads(input_file)
    #mpirun_command = ['mpirun', '-np', f"{nbeads}", os.path.join(executable_dir, 'pimdb'), '-partition', f"{nbeads}x1", '-in', input_file]
    mpirun_command = ['mpirun', '--oversubscribe', '-np', str(nbeads), (executable_dir / 'pimdb').as_posix(), '-partition', f"{nbeads}x1", '-in', input_file.as_posix()]
    
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
    
    # Check if the values of observables match
    for column in columns:
        are_equal, index = compare_arrays(data_actual[column], data_expected[column])
        if not are_equal:
            raise AssertionError(f"Test failed: '{column}' does not match at step {index}.")
        
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
            
            # Compare the output with the expected output
            if compare_output(out_path, expected_out_path):
                print("Test passed: Output matches expected output.")
            else:
                raise AssertionError("Test failed: Output does not match expected output.")
            
            # Clean up the generated output file
            print("Deleting:", out_folder)
            shutil.rmtree(out_folder)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PIMD-B++ testing of observables')
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
