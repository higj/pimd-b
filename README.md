
# PIMD-B++	


## Introduction 

Write the introduction..

## Installation

The only required library for this program is MPI (for parallelization). The code is written in accordance with the C++20 standard.

### Installing on power

To compile on HPC cluster, you need to enable the following modules:
```bash
module load cmake/cmake-3.20.0
module load mpi/openmpi-4.1.4
module load gcc/gcc-13.2.0
```
Navigate to the base directory, create a folder named build using `mkdir build`, and enter the folder. Now, configure the project and generate a native build system by running
```bash
cmake -DCMAKE_BUILD_TYPE=Release -D NDIM=1 ..
```
Note that in this example, we compile a binary for 1D systems. Change `NDIM` to the desired dimension. Currently, the following CMake cache entries are available:

* `-D NDIM=<DIM>` sets the number of spatial dimensions
* `-D OLD_BOSONIC_ALGORITHM=1` to build with the original bosonic PIMD algorithm that scales as $\mathcal{O}(PN!)$. By default, the program uses the [Feldman-Hirshberg](https://arxiv.org/abs/2305.18025) algorithm that scales as $\mathcal{O}(N^2+PN)$.

Finally, to compile the binaries, use
```bash
cmake --build .
```
## Usage

You can verify that the program is working by executing
```bash
./pimdb --dim
```
For `NDIM=1` this should print
```
Program was compiled for 1-dimensional systems.
```

### Quick Start 

To run a simulation, you will need to create a `config.ini` file. A typical config file will have the following structure:

```
; Configuration file for the PIMD simulation

[simulation]
dt = 1.0 femtosecond
steps = 1e3
threshold = 0.0
sfreq = 1
enable_thermostat = true
nbeads = 4
bosonic = true
pbc = true
seed = 12345
estimators = virial,potential
initial_position = random

[system]
temperature = 1.0 kelvin
natoms = 1
size = 1.0 picometer
mass = 1.0 atomic_unit

[interaction_potential]
name = free

[external_potential]
name = harmonic
omega = 3.0 millielectronvolt
```

Currently, the following potentials are available:

* **External potentials**: `free`, `harmonic` and `double_well` 
* **Pair interaction potentials**: `free`, `harmonic` and `dipole`

Each potential comes with its own parameters that must be provided in the configuration file. For example, the harmonic potential 
requires the `omega` parameter that sets the angular frequency of the oscillator, in energy units (e.g., `omega` is in fact $\hbar\omega$).

All pair interaction potentials require the `cutoff` parameter. By default, the cutoff distance is set to a negative value, which tells the program to calculate all interactions, regardless of the distance between a given pair of particles. Setting cutoff to zero will disable interactions altogether.
A positive cutoff aborts the calculation of the interaction force if the inter-particle distance is larger than the specified cutoff distance.

The `initial_position` option allows to specify the method of initialization for the bead coordinates. Currently, two options are available:

* `random` (default): samples random positions in each Cartesian direction from a uniform distribution on the interval $[-L/2, L/2]$, where $L$ is the linear size of the system
* `xyz(<filename>.xyz)`: initializes the coordinates based on the provided `.xyz` file. A given particle is initialized at the same location across all imaginary time-slices (beads).

Internally, the simulation uses atomic units. However, the input parameters may be provided in the units of your choosing (e.g., electron-volts for energy).

On Windows, run the program using:

```bash
mpiexec -n P pimdb.exe
```

where `P` is the number of beads.

On power, run using:
```bash
mpirun -np P pimdb -partition Px1
```

For convenience, you may use the following bash script:

```bash
#!/bin/bash

# Path to executable (by default, the name of the compiled executable is "pimdb")
EXE=/path/to/program/build/pimdb

# Path to the directory where the "output" folder will be created
RUNDIR=/path/to/simulation/folder

# Number of cores to run. In PIMD, it is usually equal to the number of beads
NCORE=P

# Clear any loaded modules
module purge

# Load modules (same modules that were used during compilation)
module load cmake/cmake-3.20.0
module load mpi/openmpi-4.1.4
module load gcc/gcc-13.2.0

# Enter the run directory
cd $RUNDIR

# Run the simulation
mpirun -np $NCORE $EXE -partition ${NCORE}x1
```

where the paths and the number of beads (`P`) must be replaced according to the specifics of your simulation. 
It can then be sent as a job using, for example,

```bash
qsub -q <queue> -lnodes=1:ppn=P -r y <script_name>.sh
```

A successful run should result in the following:

```
 __       __      __  
|__)||\/||  \ __ |__) 
|   ||  ||__/    |__)

[*] Initializing the simulation parameters
[*] Running the simulation
[*] Simulation finished running successfully
```
