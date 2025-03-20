program test_metropolis
use metropolis_monte_carlo
use statistics
implicit none

integer(kind=8) :: nstep = 10
integer :: i, nruns = 1
double precision :: energies(30), acceptances(30)
double precision :: mean_e, std_dev_e, std_dev_a, std_err_e, mean_a, std_err_a
double precision :: dt
integer :: n=1, A=1, Z(1)
double precision :: nuclei(1, 3), coeffs(1), alpha, energy, acceptance

nuclei(1,:) = [0d0, 0d0, 0d0]
!nuclei(1,:) = [0d0, 0d0, 0.7180959d0]
coeffs = 1d0
dt = 1d0
Z = 1
alpha = 1.2d0

do i = 1, nruns
    call gen_metropolis_monte_carlo(nstep, dt, n, A, Z, nuclei, coeffs, alpha, energies(i), acceptances(i))
end do
call statistical_analysis(30, energies, mean_e, std_dev_e, std_err_e)
call statistical_analysis(30, acceptances, mean_a, std_dev_a, std_err_a)

print *, "Estimated energy is: ", mean_e, " with std err: ", std_err_e
print *, "Acceptance is: ", mean_a, " with std err: ", std_err_a
end program test_metropolis
