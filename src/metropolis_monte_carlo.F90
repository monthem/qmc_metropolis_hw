module metropolis_monte_carlo
use wawefunction
use hamiltonian
use statistics
use coordinates
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

subroutine gen_metropolis_monte_carlo(nstep, dt, n, A, Z, nuclei, coeffs, alpha, energy, acceptance)
    ! Arguments
    integer(kind=8), intent(in)                         :: nstep
    double precision, intent(in)                        :: dt
    integer, intent(in)                                 :: n, A, Z(A)
    double precision, intent(in)                        :: nuclei(A, 3)
    double precision, intent(in)                        :: coeffs(A), alpha
    double precision, intent(out)                       :: energy, acceptance
    ! Local variables
    integer(kind=8)                                     :: step, accepted
    double precision                                    :: electrons(n, 3), Q_old(n, A, 3), Q_new(n, A, 3)
    double precision                                    :: chi(3), diffusion(3)
    double precision                                    :: v, a_prob, sqrt_dt, P
    double precision                                    :: drift_vectors_old(n, 3), drift_vectors_new(n, 3)
    double precision                                    :: electrons_old(n,3), electrons_new(n, 3)
    double precision                                    :: orbitals_old(n), mu_all_old(n, A),  V_NN, V_NE, V_EE, T
    double precision                                    :: orbitals_new(n), mu_all_new(n, A)
    integer                                             :: i

    ! Initialization of some variables
    energy = 0d0
    sqrt_dt = sqrt(dt)
    acceptance = 0d0
    accepted = 0_8
    ! Initialization of individual electron variables
    do i = 1, n
        electrons_old(i, :) = gaussian_rnd(3)
        electrons_new(i, :) = electrons_old(i,:) ! just to have it filled
        call form_Q(A, n, nuclei, electrons_old, Q_old)
        call MO(A, coeffs, alpha, Q_old(i,:,:), mu_all_old(i,:), orbitals_old(i))
        call calculate_drift_vector(A, coeffs, alpha, Q_old(i,:,:), mu_all_old(i,:), orbitals_old(i), drift_vectors_old(i, :))
    end do
    do step = 1, nstep
        ! START OF A MONTE CARLO RUN
        ! EACH STEP THERE ARE N ELECTRON MOVES
        do i = 1, n
           V_NE = nuc_el(A, Z, Q_old(i,:,:))
           T = kinetic_energy(A, coeffs, alpha, Q_old(i,:,:), mu_all_old(i,:), orbitals_old(i))
           energy = energy + T + V_NE
           chi = gaussian_rnd(3)
           diffusion = sqrt_dt * chi
           electrons_new(i,:) = electrons_old(i,:) + dt * drift_vectors_old(i,:) + diffusion
           call form_Q(A, n, nuclei, electrons_new, Q_new)
           call MO(A, coeffs, alpha, Q_new(i,:,:), mu_all_new(i,:), orbitals_new(i))
           call calculate_drift_vector(A, coeffs, alpha, Q_new(i,:,:), mu_all_new(i,:), orbitals_new(i), drift_vectors_new(i,:))
           P = (orbitals_new(i)/orbitals_old(i))**2
           call gmmc_acceptance_probability(electrons_old(i,:), electrons_new(i,:), drift_vectors_old(i,:), drift_vectors_new(i,:), dt, a_prob)
           a_prob = a_prob * P
           call random_number(v)
           if (v <= a_prob) then
               accepted = accepted + 1_8
               electrons_old(i,:) = electrons_new(i,:)
               Q_old(i,:,:) = Q_new(i,:,:)
               drift_vectors_old(i,:) = drift_vectors_new(i,:)
               mu_all_old(i,:) = mu_all_new(i,:)
               orbitals_old(i) = orbitals_new(i)
            end if
        end do
        V_EE = el_el(n, A, Q_new)
        energy = energy + V_EE
    end do
    energy = energy / dble(nstep)
    V_NN = nuc_nuc(A, Z, nuclei)
    energy = energy + V_NN
    acceptance = dble(accepted) / (dble(n * nstep))

end subroutine gen_metropolis_monte_carlo

end module metropolis_monte_carlo
