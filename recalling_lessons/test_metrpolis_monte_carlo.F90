program test_metropolis_monte_carlo
use section_3
implicit none

integer(kind=8) :: nmax = 100000000
double precision :: a, dt, energy, acceptance

a = 1.0d0
dt = 1.21d0

call metropolis_monte_carlo(nmax, a, dt, energy, acceptance)

print *, "Energy is: ", energy
print *, "Acceptance rate is: ", acceptance
end program test_metropolis_monte_carlo
