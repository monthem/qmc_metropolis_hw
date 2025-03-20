program test_ang2bohr
use coordinates
implicit none

integer :: A = 2
double precision :: nuclei(2, 3)

nuclei(1, :) = [1d0 , 2.5d0, 0.77d0]
nuclei(2, :) = [-1d0, -1d0, 0.74d0]

call nuclei2bohr(A, nuclei)

print *, nuclei
end program test_ang2bohr
