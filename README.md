# PIMD-B++	
![Ubuntu compilation check](https://github.com/higj/pimd-b/actions/workflows/cmake-ubuntu.yml/badge.svg)

## Introduction 

### Bosonic algorithm

#### Old algorithm

In the case of bosons, the partition function must be evaluated in a properly symmetrized basis. As a result, the partition function contains contributions from $N!$ permutations, where $N$ is the number of quantum particles. A naive bosonic PIMD algorithm takes into account all the $N!$ permutations. Labeling the spring energy due to a specific permutation $\sigma$ as

$$
E^{\sigma}=\frac{1}{2}m\omega_{P}^{2}\sum_{\ell=1}^{N}\sum_{j=1}^{P}\left(\mathbf{r}_ {\ell}^{j}-\mathbf{r}_ {\ell}^{j+1}\right)^{2},\\;\text{where} \\; \mathbf{r}_ {\ell}^{P+1}=\mathbf{r}_ {\sigma\left(\ell\right)}^{1},
$$

one can define an effective bosonic ring polymer potential as

$$
V_B =-\frac{1}{\beta}\ln\left[\frac{1}{N!}\sum_{\sigma}e^{-\beta E^{\sigma}}\right].
$$

The spring force on a given bead $\mathbf{r}_{\ell}^{j}$ is then given by

$$
\mathbf{f}_ {\ell}^{j}=-\nabla_{\mathbf{r}_ {\ell}^{j}}V_{B}=-\frac{\sum_{\sigma}e^{-\beta E^{\sigma}}\nabla_{\mathbf{r}_ {\ell}^{j}}E^{\sigma}}{\sum_{\sigma}e^{-\beta E^{\sigma}}}.
$$

It is not actually necessary to evaluate all the terms in $E^{\sigma}$ when calculating the weights. Indeed, let us define $\Delta E^{\sigma}$ as 

$$
\Delta E^{\sigma}=E^{\sigma}-\frac{1}{2}m\omega_ {P}^{2}\sum_ {\ell=1}^{N}\sum_ {j=1}^{P-1}\left(\mathbf{r}_ {\ell}^{j}-\mathbf{r}_ {\ell}^{j+1}\right)^{2}=\frac{1}{2}m\omega_ {P}^{2}\sum_{\ell=1}^{N}\left(\mathbf{r}_ {\ell}^{P}-\mathbf{r}_ {\sigma\left(\ell\right)}^{1}\right)^{2}.
$$

Since $E^{\sigma} - \Delta E^{\sigma}$ is the same for all permutations, we can write

$$
V_{B}=\left(E^{\sigma}-\Delta E^{\sigma}\right)-\frac{1}{\beta}\ln\left[\frac{1}{N!}\sum_{\sigma}e^{-\beta\Delta E^{\sigma}}\right],
$$

and also

$$
\mathbf{f}_ {\ell}^{j}=-\frac{\sum_{\sigma}e^{-\beta\Delta E^{\sigma}}\nabla_{\mathbf{r}_ {\ell}^{j}}E^{\sigma}}{\sum_{\sigma}e^{-\beta\Delta E^{\sigma}}}.
$$

For interior beads, that is, for $j=2,\dots,P-1$, the forces coincide with the ordinary spring forces in systems of distinguishable particles, i.e.,

$$
-\nabla_{\mathbf{r}_ {\ell}^{j}}E^{\sigma}=-m\omega_ {P}^{2}\left(2\mathbf{r}_ {\ell}^{j}-\mathbf{r}_ {\ell}^{j+1}-\mathbf{r}_ {\ell}^{j-1}\right).
$$

Bosonic exchange affects only the forces acting on the exterior beads ($j=1,P$). Note however that the aforementioned optimizations do not change the scaling of this algorithm which remains $\mathcal{O}(N!)$.

For the primitive kinetic energy estimator, we use the following convenient expression:

$$
\left\langle K\right\rangle =\frac{dNP}{2\beta}-E_{\mathrm{sp,interior}}-\left\langle \frac{\sum_{\sigma}\Delta E^{\sigma}e^{-\beta\Delta E^{\sigma}}}{\sum_{\sigma}e^{-\beta\Delta E^{\sigma}}}\right\rangle,
$$

where

$$
E_{\mathrm{sp,interior}}=\frac{1}{2}m\omega_{P}^{2}\sum_{\ell=1}^{N}\sum_{j=1}^{P-1}\left(\mathbf{r}_ {\ell}^{j}-\mathbf{r}_ {\ell}^{j+1}\right)^{2},
$$

is the spring energy of the interior beads.


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
* **Pair interaction potentials**: `free`, `harmonic`, `dipole` and `aziz`

Each potential comes with its own parameters that must be provided in the configuration file. For example, the harmonic potential 
requires the `omega` parameter that sets the angular frequency of the oscillator, in energy units (i.e., `omega` is in fact $\hbar\omega$).

All pair interaction potentials require the `cutoff` parameter. By default, the cutoff distance is set to a negative value, which tells the program to calculate all interactions, regardless of the distance between a given pair of particles. Setting cutoff to zero will disable interactions altogether.
A positive cutoff aborts the calculation of the interaction force if the inter-particle distance is larger than the specified cutoff distance.

The `initial_position` option allows to specify the method of initialization for the bead coordinates. Currently, two options are available:

* `random` (default): samples random positions in each Cartesian direction from a uniform distribution on the interval $[-L/2, L/2]$, where $L$ is the linear size of the system
* `xyz(<filename>.xyz)`: initializes the coordinates based on the provided `.xyz` file. A given particle is initialized at the same location across all imaginary time-slices (beads).

The `size` option defines the linear size of the system. Currently, only cube geometry is supported. In the absence of periodic boundary conditions, `size` only affects the way initial positions are generated. However, if periodic boundary conditions are enabled, the system size also affects the cutoff distance for interactions, as well as the estimators. Also, the coordinates may be wrapped in this case, and minimum image convention can potentially be employed, if such functionality is desired.

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
|`seed`     | Random number generator seed |
|`initial_position`     | Method for generating the initial positions of the beads |
|`initial_velocity`     | Method for generating the initial velocities of the beads. `random` samples from the Maxwell-Boltzmann distribution. `manual` loads velocities from a provided file. (Default: `random`) |
|`temperature`     |  Temperature of the system (units of temperature) |
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
mpirun -np P pimdb -partition Px1
```

If the configuration file for the simulation is not located in the same directory as the executable, you can specify the path to the configuration file using the `-in` flag, e.g.,

```bash
mpirun -np P pimdb -partition Px1 -in /path/to/config.ini
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
