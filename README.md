## Description  
This program computes the electronic energies of small **homonuclear systems** (e.g., H, He, H₂⁺, H₃⁺, Li⁺, Be²⁺) using a Variational Monte Carlo (VMC) method with a generalized Metropolis algorithm. The wavefunction ansatz is a Hartree product of molecular orbitals (MOs) formed as Linear Combinations of Atomic Orbitals (LCAO) with 1s Slater basis functions. The user specifies the system geometry, basis set exponents, LCAO coefficients, and simulation parameters in an input file.

## Directory Structure  
- **`/src`**: Contains all Fortran source code.  
  - `file_IO.F90`: Reads and parses input files.  
  - `coordinates.F90`: Manages nuclear and electron coordinates.  
  - `statistics.F90`: Computes statistical values and generates Gaussian-distributed random numbers.  
  - `wawefunction.F90`: Handles wavefunction, ∇ψ, and ∇²ψ calculations.  
  - `hamiltonian.F90`: Evaluates Hamiltonian components.  
  - `metropolis_monte_carlo.F90`: Implements the generalized Metropolis algorithm.  
  - `main.F90`: Main driver for the simulation.  
- **`/examples`**: Example input files for testing.

## Compilation  
1. Navigate to the `src` directory:  
   ```bash  
   cd src
2. Compile the program:
    ```bash
    gfortran file_IO.F90 coordinates.F90 statistics.F90 wawefunction.F90 hamiltonian.F90 metropolis_monte_carlo.F90 main.F90 -o monte_carlo

## Usage
1. Run the program:
    ```bash
    /path/to/executable/ /path/to/input.inp
    ./executable input.inp
    ./executable input.inp > output.out
## Example Output
    ```bash
    $ ./monte_carlo example.inp
Input file read succesfully!
Performing           30  Monte Carlo runs now.
Monte Carlo runs finished succesfully!
Computation time:                   00min 38sec
Computed electronic energy:         -0.435991689 a.u.
Standard deviation:                     0.001386 a.u.
Standard error:                         0.000253 a.u.

Acceptance rate:                          0.6208
Standard deviation:                       0.0005
Standard error:                           0.0001

## Input File Format
SYSTEM
<N_atoms> <N_electrons> <N_atomtypes>
<Z> <x> <y> <z> [repeated for each atom]
BASIS
<Z> <alpha> [repeated for each atom type]
LCAO
<c1> <c2> ... <cN_atoms>
SIMULATION
<N_runs> <N_steps> <dt>
COMMENT
[Optional user notes]

## Requirements
GNU Fortran compiler (gfortran). Program was tested with gfortran 14.2.1.

For installation instructions, refer to INSTALL.md
