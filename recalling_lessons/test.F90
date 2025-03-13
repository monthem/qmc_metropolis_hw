program test
implicit none

double precision :: r1(3), r2(3)

r1 = [1d0, 2d0, 3d0]
r2 = [2d0, 2d0, 2d0]

print *, r1 * r2
print *, exp(r1)
end program test
