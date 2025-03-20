# Installation Instructions  

Installation and usage instructions assume Linux based OS.

## Requirements  
- **GNU Fortran compiler** (`gfortran`).  

## Installing `gfortran`  

### Ubuntu/Debian:  
    sudo apt install gfortran 
### Arch Linux:
    sudo pacman -S gcc-fortran

## Compilation
1. Navigate to the `src` directory:
   ```bash
   cd src
2. Compile the program:
    ```bash
    gfortran file_IO.F90 coordinates.F90 statistics.F90 wawefunction.F90 hamiltonian.F90 metropolis_monte_carlo.F90 main.F90 -o monte_carlo

