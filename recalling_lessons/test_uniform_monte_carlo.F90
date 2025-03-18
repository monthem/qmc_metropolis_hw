program test_uniform_monte_carlo
use section_3
implicit none

double precision        :: a = 1.0d0
integer(kind=8)         :: nmax = 100000000 
double precision        :: energy 

energy = uniform_monte_carlo(nmax, a)

print *, energy

end program test_uniform_monte_carlo
