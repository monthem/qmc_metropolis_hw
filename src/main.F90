program main_program
!use omp_lib
use file_IO
use coordinates
use statistics
use wawefunction
use hamiltonian
use metropolis_monte_carlo
implicit none

! Input file variables
character(len=256)                          :: filename
integer                                     :: n, A, N_A, N_runs
integer, allocatable                        :: Z(:)
double precision, allocatable               :: R(:,:), alphas(:), coeffs(:)
double precision                            :: dt, alpha
integer(kind=8)                             :: N_step
! Statistics variables
double precision, allocatable               :: energies(:), acceptances(:)
double precision                            :: mean_e, mean_a, sdev_e, sdev_a, serr_e, serr_a
! Local variables
integer                                     :: i, minutes, seconds
double precision                            :: t1, t2, runtime

! Check if path to an input file was provided
if (command_argument_count() /= 1) then
    print *, "Usage: ./program_name inputfile.inp"
    stop
end if
! Get the provided input file name
call get_command_argument(1, filename)

! Read the input file
call read_input_file(filename, A, n, N_A, Z, R, alphas, coeffs, N_runs, N_step, dt)
! Get the alpha value
alpha = alphas(1)                           ! Hetero diatomics not supported yet. For now only one alpha value.
! Convert the read geometry from Angstrom to Bohr
call nuclei2bohr(A, R)

allocate(energies(N_runs), acceptances(N_runs))

print *, "Input file read succesfully!"
print *, "Performing ", N_runs, " Monte Carlo runs now."
! Perform N_runs Monte Carlo runs
call cpu_time(t1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
do i = 1, N_runs
    call gen_metropolis_monte_carlo(N_step, dt, n, A, Z, R, coeffs, alpha, energies(i), acceptances(i))
end do
!$OMP END PARALLEL DO
call cpu_time(t2)

! Calculate runtime
runtime = t2 - t1
minutes = int(runtime / 60d0)
seconds = mod(int(runtime), 60)

! Get the statistics for the N_runs
call statistical_analysis(N_runs, energies, mean_e, sdev_e, serr_e)
call statistical_analysis(N_runs, acceptances, mean_a, sdev_a, serr_a)

print *, "Monte Carlo runs finished succesfully!"
print '(A, I2.2, A, I2.2, A)', "Computation time:                   ", minutes, "min ", seconds, "sec"

! Print the output values
print '(A, F18.9, A)', "Computed electronic energy:   ", mean_e, " a.u."
print '(A, F18.6, A)', "Standard deviation:           ", sdev_e, " a.u."
print '(A, F18.6, A)', "Standard error:               ", serr_e, " a.u."
print *

print '(A, F18.4)', "Acceptance rate:              ", mean_a
print '(A, F18.4)', "Standard deviation:           ", sdev_a
print '(A, F18.4)', "Standard error:               ", serr_a

deallocate(Z, R, alphas, coeffs, energies, acceptances)

end program main_program
