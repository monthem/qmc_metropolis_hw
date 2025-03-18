program test_gradient_vector
use wawefunction
implicit none

integer :: N = 2
double precision :: coeffs(2), alpha, q(2,3), mu(2)
double precision :: phi
double precision :: gradient(3)



print *, "testing grad vector now"

print *, "H atom with c = 1.0 and alpha = 1.2"

N = 2
coeffs = [1.0d0, 1.0d0]
alpha = 1.2d0
q(1, :) = [-1.d0, 0.5d0, -0.8d0]
q(2, :) = [-1d0, 1d0, 1d0]
call MO(N, coeffs, alpha, q, mu, phi)

print *, "H atom: MO (ie the WF)"
print *, "MO is: ", phi

call gradient_vector(N, coeffs, alpha, q, mu, gradient)

print *, "gradient vector is: ", gradient
print *, "drift vector is: ", gradient / phi
end program test_gradient_vector
