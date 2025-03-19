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
