module file_IO
implicit none
contains

subroutine read_input_file(filename, A, n, N_A, Z, R, alphas, coeffs, N_runs, N_step, dt)
    ! Arguments
    character(len=*), intent(in)                :: filename
    integer, intent(out)                        :: A, N_A, n, N_runs
    integer, allocatable, intent(out)           :: Z(:)
    double precision, allocatable, intent(out)  :: R(:,:), alphas(:), coeffs(:)
    double precision, intent(out)               :: dt
    integer(kind=8), intent(out)                :: N_step
    ! Local variables
    integer                                     :: i, j, ios, unit
    character(len=20)                           :: section

    ! Open the input file
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) stop "Error opening input file. Confirm that file exists or is written properly."
    ! Read SYSTEM section
    read(unit, *) section
    if (trim(section) /= "SYSTEM") stop "Expected SYSTEM section. Input file is incorrect!"
    read(unit, *) A, n, N_A
    allocate(Z(A), R(A, 3), alphas(N_A), coeffs(A))
    do i = 1, A
        read(unit, *, iostat=ios) Z(i), R(i, :)
    end do
    ! Read BASIS section
    read(unit, *) section
    if (trim(section) /= "BASIS") stop "Expected BASIS section. Input file is incorrect!"
    do i = 1, N_A
        read(unit, *) j, alphas(i)
    end do
    ! Read LCAO section
    read(unit, *) section
    if (trim(section) /= "LCAO") stop "Expected LCAO section. Input file is incorrect!"
    read(unit, *) coeffs(:)
    ! Read SIMULATION section
    read(unit, *) section
    if (trim(section) /= "SIMULATION") stop "Expected SIMULATION. Input file is incorrect!"
    read(unit, *) N_runs, N_step, dt
    ! Close file
    close(unit)
end subroutine read_input_file

end module file_IO
