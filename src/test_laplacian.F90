program test_laplacian
use wawefunction
implicit none

integer :: N = 1
double precision :: coeffs(1), alpha, q(1,3), mu(1)
double precision :: phi
double precision :: del_sq



print *, "testing Laplacian"

print *, "H atom with c = 1.0 and alpha = 1.2"

N = 1
coeffs = [1.0d0]
alpha = 1.2d0
q(1, :) = [-1.d0, 0.5d0, -0.8d0]
!q(2, :) = [-1d0, 1d0, 1d0]
call MO(N, coeffs, alpha, q, mu, phi)

!print *, "H atom: MO (ie the WF)"
!print *, "MO is: ", phi

call laplacian(N, coeffs, alpha, q, mu, del_sq)

print *, "Laplacian is is: ", del_sq
print *, "kinetic energy is: ", -0.5d0 * del_sq / phi
end program test_laplacian
