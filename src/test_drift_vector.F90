program test_drift_vector
use wawefunction
use metropolis_monte_carlo
implicit none

integer :: A
double precision :: coeffs(1), alpha = 1.2d0, mu(1), phi, drift_vector(3)
double precision :: q(1, 3)

A = 1
coeffs = [1d0]
q(1, :) = [-1d0, 0.5d0, -0.8d0]
!q(2, :) = [-0.55d0, 0d0, 0.25d0]
call MO(A, coeffs, alpha, q, mu, phi)
call calculate_drift_vector(A, coeffs, alpha, q, mu, phi, drift_vector)

print *, drift_vector

end program test_drift_vector
