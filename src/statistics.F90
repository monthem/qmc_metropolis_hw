module statistics
implicit none
contains

subroutine statistical_analysis(N, values, mean, std_dev, std_err)
    ! Arguments
    integer, intent(in)              :: N
    double precision, intent(in)     :: values(N)
    double precision, intent(out)    :: mean, std_dev, std_err
    ! Local variable
    double precision                 :: variance

    mean = sum(values) / dble(N)
    variance = sum((values-mean)**2) / dble(N)
    std_dev = sqrt(variance)
    std_err = std_dev / sqrt(dble(N))
end subroutine statistical_analysis

subroutine generate_rnd_r(r, a, b)
    ! Generates a vector r with components being uniform random numbers in range [a,b]

    ! Arguments
    double precision, intent(in)        :: a, b
    double precision, intent(out)       :: r(3)
    ! Local variables
    double precision                    :: x, y, z

    ! Generate three random numbers between 0 and 1
    call random_number(x)
    call random_number(y)
    call random_number(z)

    ! Scale them to be in wanted range such that if rand x is 1, x will be max[a,b] and if rand x is 0, x is min[a,b]
    x = x * (b - a) + a
    y = y * (b - a) + a
    z = z * (b - a) + a

    r = [x, y, z]
end subroutine generate_rnd_r

function gaussian_rnd(n) result(z)
    ! Argument
    integer, intent(in) :: n
    ! Output
    double precision :: z(n)
    ! Local variable
    double precision :: u(n+1)
    double precision, parameter :: two_pi = 2.0d0*acos(-1.d0)
    integer :: i

    call random_number(u)
    if (iand(n, 1) == 0) then
        do i = 1, n, 2
            z(i) = sqrt(-2.d0*log(u(i)))
            z(i+1) = z(i) * sin(two_pi*u(i+1))
            z(i) = z(i) * cos(two_pi*u(i+1))
        end do
    else
        do i = 1, n-1, 2
            z(i) = sqrt(-2.d0*log(u(i)))
            z(i+1) = z(i) * sin(two_pi*u(i+1))
            z(i) = z(i) * cos(two_pi*u(i+1))
        end do

        z(n) = sqrt(-2.d0*log(u(n)))
        z(n) = z(n) * cos(two_pi*u(n+1))
    end if
end function gaussian_rnd


end module statistics
