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
        call generate_rnd_r(u, -1d0, 1d0) ! NO! CORRECT TOMORROW!
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

subroutine calculate_drift_vector(a, r, d)
    ! drift vector for H atom wawefunction
    ! Arguments
    double precision, intent(in)    :: a, r(3)
    double precision, intent(out)   :: d(3)
    ! Local variables
    double precision :: psi, del_psi(3), temp_var
    psi = wawefunction(a, r)
    temp_var = -1.0d0 * a / sqrt(sum(r**2))
    del_psi = [r(1)*temp_var, r(2)*temp_var, r(3)*temp_var]

    d = del_psi / psi
end subroutine calculate_drift_vector

subroutine gmmc_acceptance_probability(r_old, r_new, d_old, d_new, dt, a_prob)
    ! Arguments
    double precision, intent(in)        :: r_old(3), r_new(3)
    double precision, intent(in)        :: d_old(3), d_new(3)
    double precision, intent(in)        :: dt
    double precision, intent(out)       :: a_prob

    a_prob = (r_new - r_old) * (0.5d0 * dt * (d_new**2 - d_old**2))
    a_prob = exp(a_prob)    ! Do not forget to multiply by |psi(r')|² / |psi(r)|² in gmmc subroutine!
end subroutine gmmc_acceptance_probability

subroutine gen_metropolis_monte_carlo(a, nmax, dt, energy, acceptance)
    ! Arguments
    integer(kind=8)                         :: nmax
    double precision, intent(in)            :: a, dt
    double precision, intent(out)           :: energy, acceptance
    ! Local variables
    integer(kind=8)                         :: n, accepted = 0
    double precision                        :: r_old(3), r_new(3), psi_old, psi_new
    double precision                        :: d_old(3), d_new(3), chi(3)
    double precision                        :: v, a_prob, nmax_inv, diffusion, sqrt_dt, P
    
    ! Initialization of some variables
    acceptance = 0d0
    u = 0d0
    nmax_inv = 1d0 / dble(nmax)
    sqrt_dt = sqrt(dt)
    call generate_rnd_r(r_old, -5d0, 5d0)
    energy = local_energy(a, r_old)

    do n = 2, nmax
        chi = gaussian_rnd(3)
        diffusion = sqrt_dt * chi
        call calculate_drift_vector(a, r_old, d_old)
        r_new = r_old + dt * d_old + diffusion
        call calculate_drift_vector(a, r_new, d_new)
        energy = energy + local_energy(a, r_new)
        psi_old = wawefunction(a, r_old)
        psi_new = wawefunction(a, r_new)
        P = (psi_new / psi_old)**2
        call gmmc_acceptance_probability(r_old, r_new, d_old, d_new, dt, a_prob)
        a_prob = a_prob * P
        call random_number(v)
        if (v <= a_prob) then
            r_old = r_new
            accepted = accepted + 1
        end if
    end do
    energy = energy * nmax_inv
    acceptance = accepted * nmax_inv
end subroutine gen_metropolis_monte_carlo

end module section_3
