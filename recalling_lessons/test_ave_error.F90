program test_ave_error
use section_3
implicit none

integer         :: n = 3
double precision :: x(3)
double precision :: ave, error

x = [1, 2, 3]

call ave_error(n, x, ave, error)
print *, ave, error
end program test_ave_error
