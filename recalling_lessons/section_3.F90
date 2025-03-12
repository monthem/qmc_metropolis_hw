module section_3
use H_local_energy
implicit none
contains

subroutine ave_error(n, x, ave, error)
    ! Arguments
    integer                         :: n
    double precision, intent(in)    :: x(n)
    double precision, intent(out)   :: ave, error
    ! Local variable
    double precision                :: variance

    if (n < 1) then
        stop "Error: Empty array in ave_error subroutine!)"
    else if (n==1) then
        ave = x(1)
        error = 0.0d0
    else
        ave = sum(x) / dble(n)
        variance = sum((x - ave)**2) / dble(n-1)
        error = sqrt(variance/dble(n))
    end if
end subroutine ave_error

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


function uniform_monte_carlo(nmax, a) result(energy)
    ! Arguments
    integer(kind=8), intent(in)         :: nmax
    double precision, intent(in)        :: a
    ! Output
    double precision                    :: energy
    ! Local variables
    double precision                    :: normalization, r(3)
    integer(kind=8)                     :: n

    normalization = 0d0
    energy = 0d0
    
    do n = 1, nmax
        r = 0.0d0
        call generate_rnd_r(r, -5d0, 5d0)
        normalization = normalization + wawefunction(a, r)**2
        energy = energy + local_energy(a, r) * wawefunction(a, r)**2
    end do

    energy = energy / normalization
end function uniform_monte_carlo

subroutine metropolis_monte_carlo(nmax, a, dt, energy, acceptance)
    ! Performs a monte carlo integration with metropolis sampling
    
    ! Arguments
    integer(kind=8)                 :: nmax
    double precision, intent(in)    :: a, dt
    double precision, intent(out)   :: energy, acceptance
    ! Local variables
    integer(kind=8)                 :: n, accepted = 0
    double precision                :: r_old(3), r_new(3), u(3), psi_old, psi_new
    double precision                :: ratio, v, nmax_inv

    ! Initialization of some variables
    acceptance = 0d0
    call generate_rnd_r(r_old, -5d0, 5d0)
    u = 0d0
    nmax_inv = 1d0 / dble(nmax)
    energy = local_energy(a, r_old)

    do n = 2, nmax
        call generate_rnd_r(u, -1d0, 1d0)
        r_new = r_old + dt * u
        energy = energy + local_energy(a, r_new)
        psi_old = wawefunction(a, r_old)
        psi_new = wawefunction(a, r_new)
        ratio = (psi_new / psi_old)**2
        call random_number(v)
        if (v <= ratio) then
            r_old = r_new
            accepted = accepted + 1
        end if
    end do
    energy = energy * nmax_inv
    acceptance = accepted * nmax_inv
end subroutine metropolis_monte_carlo
        


    
    

end module section_3
