program test_inp
use file_IO
implicit none

character(len=256) :: filename
integer :: A, n, N_A, N_runs
integer, allocatable :: Z(:)
double precision, allocatable :: R(:, :), alphas(:), coeffs(:)
double precision :: dt
integer(kind=8) :: N_step

integer :: i

filename = "example_input.inp"

call read_input_file(filename, A, n, N_A, Z, R, alphas, coeffs, N_runs, N_step, dt)

print *, A, N, N_A
do i = 1, A
    print *, Z(i), R(A, :), "alpha for this nucleus: ", alphas
end do
print *, coeffs
print *, N_runs, N_step, dt

end program test_inp

