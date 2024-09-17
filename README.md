# PIMD-B++	
![Ubuntu compilation check](https://github.com/higj/pimd-b/actions/workflows/cmake-ubuntu.yml/badge.svg)

## Introduction 

**PIMD-B++** is a lightweight program for performing Boltzmannonic and bosonic path integral molecular dynamics simulations. For the theory behind the bosonic algorithm consult the [Wiki](https://github.com/higj/pimd-b/wiki/Bosonic-algorithm) page.

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
* `-D OLD_BOSONIC_ALGORITHM=1` to build with the original bosonic PIMD algorithm that scales as $\mathcal{O}(N!)$. By default, the program uses the [Feldman-Hirshberg](https://arxiv.org/abs/2305.18025) algorithm that scales as $\mathcal{O}(N^2+PN)$.

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

```ini
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

[output]
positions = false
velocities = false
forces = false

[observables]
energy = kelvin
classical = off
bosonic = off
```

Currently, the following potentials are available:

* **External potentials**: `free`, `harmonic` and `double_well` 
* **Pair interaction potentials**: `free`, `harmonic`, `dipole` and `aziz`

Each potential comes with its own parameters that must be provided in the configuration file. For example, the harmonic potential 
requires the `omega` parameter that sets the angular frequency of the oscillator, in energy units (i.e., `omega` is in fact $\hbar\omega$).

All pair interaction potentials require the `cutoff` parameter. By default, the cutoff distance is set to a negative value, which tells the program to calculate all interactions, regardless of the distance between a given pair of particles. Setting cutoff to zero will disable interactions altogether.
A positive cutoff aborts the calculation of the interaction force if the inter-particle distance is larger than the specified cutoff distance.

The `initial_position` option allows to specify the method of initialization for the bead coordinates. Currently, three options are available:

* `random` (default): samples random positions in each Cartesian direction from a uniform distribution on the interval $[-L/2, L/2]$, where $L$ is the linear size of the system.
* `xyz(<filename>.xyz)`: initializes the coordinates based on the provided `.xyz` file. A given particle is initialized at the same location across all imaginary time-slices (beads).
* `xyz(<filename_format>)`: if the provided filename is a [Python format string](https://docs.python.org/3/library/string.html#formatspec), the indices of the imaginary time-slices (starting with either 0 or 1, automatically detected from the available files) are substituted as the format argument, and the resulting filenames are then used to initialize the coordinates. The formatted string can contain only a single replacement field.

Similarly, the `initial_velocity` option gives the user the ability to initialize the bead velocities. Three options are available as of now:

* `random` (default): samples velocities from the Maxwell-Boltzmann distribution at the given temperature of the simulation.
* `manual`: the velocities are initialized from the files `init/vel_XX.dat` where `XX` symbolizes the two-digit index of the imaginary time-slice, starting from 01.
* `manual(<filename_format>)`: similar behavior to `xyz(<filename_format>)`.

The `size` option defines the linear size of the system. Currently, only cube geometry is supported. In the absence of periodic boundary conditions, `size` only affects the way initial positions are generated. However, if periodic boundary conditions are enabled, the system size also affects the cutoff distance for interactions, as well as the estimators. Also, the coordinates may be wrapped in this case, and minimum image convention can potentially be employed, if such functionality is desired.

In the `[output]` section, users can request the output of various quantities related to the state of the system (such as positions, velocities, etc.) 
The format for this section is `state_name = state_unit`. The key (`state_name`) must correspond to a name of a supported state. 
The value (`state_unit`) specifies the unit to which the output must be converted. If set to `false` (or, equivalently, `off`), the state 
will not be printed. By default, all states are set to `false`. If set to `true` (or, equivalently, `on`), the state will be printed in default (atomic) units, assuming the quantity is not dimensionless. Otherwise, the user 
must specify the desired unit.

In the `[observables]` section, users can specify which physical observables should be evaluated and printed in `simulation.out`. The format for this section is `observable_name = observable_unit`. The key (`observable_name`) must correspond to a name of a supported observable. For the value (`observable_unit`), users can indicate the units to which the results should be converted (if the observable is not dimensionless), or use `off` if the observable should not be calculated at all (this is the default setting for all observables except `energy`). For dimensionless estimators, users may leave the value empty or specify `none` as the unit.

Currently, three observable *types* are supported:

* `energy`: Calculates the quantum energy of the system using different estimators. Currently, the thermodynamic (primitive), virial, and potential energy estimators are supported.
* `classical`: Calculates observables related to the classical ring-polymer system, such as the kinetic energy (due to the fictitious momenta), spring energies, and temperature.
* `bosonic`: Calculates the probabilities of two types of topologies: where all particles are separate and where all particles are connected (dimensionless estimator). Printed *only* in bosonic simulations.

Internally, the simulation uses atomic units. However, the input parameters may be provided in the units of your choosing (e.g., electron-volts for energy).

The following options are available in the `[simulation]`, `[system]` and `[output]` sections of the configuration file:

| Option | Description |
| :-----------: | ----------- |
|`dt`     |  Timestep for MD (units of time) |
|`steps`     |  Total number of MD steps (can be written in scientific notation) |
|`threshold`     |  Defines the percentage of steps to throw away for thermalization (float between 0 and 1) |
|`sfreq`     | Frequency at which the observables are calculated in the production stage |
|`enable_thermostat`     |  If set to `false`, the Langevin thermostat will be disabled (NVE run) |
|`gamma`    |  Friction coefficient for the Langevin thermostat in units of inverse time (Default: $\frac{1}{100\Delta t}$) |
|`nbeads`     |  Number of imaginary time-slices (beads) |
|`bosonic`     |  Set to `true`/`false` for bosonic/distinguishable PIMD (Default: `false`) |
|`pbc`     |  Set to `true` to enable periodic boundary conditions (Default: `false`) |
|`fixcom`     |  Set to `true` to remove the center of mass motion (Default: `true`) |
|`seed`     | Random number generator seed (a positive integer below $9 \times 10^8$) |
|`initial_position`     | Method for generating the initial positions of the beads |
|`initial_velocity`     | Method for generating the initial velocities of the beads. `random` samples from the Maxwell-Boltzmann distribution. `manual` loads velocities from a provided file. (Default: `random`) |
|`temperature`     |  Temperature of the quantum system (units of temperature) |
|`natoms`     |  Number of particles in the quantum system |
|`size`     |  Linear size of the system (units of length) |
|`mass`     |  Mass of the particles (units of mass) |
|`positions`     | Set to `true` to output the trajectories (Default: `false`) |
|`velocities`     | Set to `true` to output the velocities (Default: `false`) |
|`forces`     | Set to `true` to output the forces (Default: `false`) |

On Windows, run the program using:

```bash
mpiexec -n P pimdb.exe
```

where `P` is the number of beads.

On power, run using:
```bash
mpirun -np P pimdb
```

If the configuration file for the simulation is not located in the same directory as the executable, you can specify the path to the configuration file using the `-in` flag, e.g.,

```bash
mpirun -np P pimdb -in /path/to/config.ini
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
mpirun -np $NCORE $EXE
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