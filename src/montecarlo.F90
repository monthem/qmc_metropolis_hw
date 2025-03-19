module metropolis_monte_carlo
use wawefunction
use hamiltonian
use statistics
implicit none
contains

subroutine calculate_drift_vector(A, coeffs, alpha, q, mu, phi, drift_vector)
    ! Arguments 
    integer, intent(in)                 :: A
    double precision, intent(in)        :: coeffs(A), alpha, q(A, 3), mu(A), phi
    double precision, intent(out)       :: drift_vector(3)
    ! Local variables
    double precision                    :: gradient(3)

    call gradient_vector(A, coeffs, alpha, q, mu, gradient)
    drift_vector = gradient / phi
end subroutine calculate_drift_vector

end module metropolis_monte_carlo

subroutine gmmc_acceptance_probability(r_old, r_new, d_old, d_new, dt, a_prob)
    ! Arguments
    double precision, intent(in)        :: r_old(3), r_new(3)
    double precision, intent(in)        :: d_old(3), d_new(3)
    double precision, intent(in)        :: dt
    double precision, intent(out)       :: a_prob

    a_prob = dot_product((r_new - r_old), (d_new + d_old))
    a_prob = a_prob + 0.5d0 * dt * (sum(d_new**2) - sum(d_old**2))
    a_prob = exp(-1.0d0 * a_prob)
end subroutine gmmc_acceptance_probability

