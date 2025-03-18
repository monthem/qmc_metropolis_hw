program test_gen_metropolis
use section_3
implicit none

integer(kind=8) :: nmax = 1000000
double precision :: a, dt, energy, acceptance

a = 1.2d0
dt = 1.0d0

call gen_metropolis_monte_carlo(a, nmax, dt, energy, acceptance)

print *, "Energy is: ", energy
print *, "Acceptance rate is: ", acceptance
end program test_gen_metropolis
