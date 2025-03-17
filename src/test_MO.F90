program test_AO
use wawefunction
implicit none

integer :: N = 2
double precision :: coeffs(2), alphas(2), q(2,3), mu(2)
double precision :: phi

coeffs = [1.0d0, 1.0d0]
alphas = [1.0d0, 1.0d0]

q(1, :) = [-0.5d0, 0.0d0, 0.5d0]
q(2, :) = [-1.22d0, 1.22d0, 0.0d0]

mu = 0d0

call AO(N, coeffs, alphas, q, mu)

print *, mu

call MO(N, coeffs, alphas, q, mu, phi)

print *, phi

end program test_AO
