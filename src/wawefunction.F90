module wawefunction     ! functions for the calculation of individual AOs mu, MOs phi and the total WF psi
implicit none
contains

subroutine AO(A, coeffs, alpha, q, mu)
    ! subroutine that takes the number of nuclei A, an array of AO coefficients and exponents of size N, an array of position vectors q of size N by 3, and returns an array of values of the AOs mu of size N
    ! Arguments
    integer, intent(in)             :: A 
    double precision, intent(in)    :: coeffs(A), alpha, q(A, 3)
    double precision, intent(out)   :: mu(A)
    ! Local variable
    integer                         :: i

    do i = 1, A
        mu(i) = coeffs(i) * exp(-alpha*sqrt(sum(q(i,:)**2)))
    end do
end subroutine AO

subroutine MO(A, coeffs, alpha, q, mu, phi)
    ! Arguments
    integer, intent(in)                 :: A
    double precision, intent(in)        :: coeffs(A), alpha, q(A, 3)
    double precision, intent(inout)     :: mu(A)
    double precision, intent(out)       :: phi
    ! Local variable
    integer                             :: i

    call AO(A, coeffs, alpha, q, mu)
    phi = sum(mu)
end subroutine MO

function WF(MO_1, MO_2) result(psi)
    ! Arguments
    double precision, intent(in)        :: MO_1, MO_2
    ! Output 
    double precision                    :: psi
    psi = MO_1 * MO_2
end function WF


subroutine gradient_vector(A, coeffs, alpha, q, mu, gradient)
    ! WARNING: Correct gradient vector is additionally multiplied by other (if they exist) electrons phi values
    !          Here it is going to be used only for the drift vector, which divides the gradient vector by Psi
    !          which cancels other electron's phi values
    ! Arguments
    integer, intent(in)                 :: A
    double precision, intent(in)        :: coeffs(A), alpha, q(A, 3), mu(A)
    double precision, intent(out)       :: gradient(3)
    ! Local variables
     integer                             :: i
    double precision                    :: norm_inv
    
    gradient = 0d0
    do i = 1, A
        norm_inv = 1d0 / sqrt(sum(q(i, :)**2))
        gradient = gradient + coeffs(i)*mu(i)*norm_inv*q(i,:)
    end do
    gradient = -1d0 * alpha * gradient
end subroutine gradient_vector

subroutine laplacian(A, coeffs, alpha, q, mu, del_sq)
    ! Arguments
    integer, intent(in)                 :: A
    double precision, intent(in)        :: coeffs(A), alpha, q(A, 3), mu(A)
    double precision, intent(out)       :: del_sq
    ! Local variables
    integer                             :: i
    double precision                    :: norm_inv
    
    del_sq = 0d0
    do i = 1, A
        norm_inv = 1d0 / sqrt(sum(q(i, :)**2))
        del_sq = del_sq + coeffs(i)*(alpha**2 - 2d0 * alpha * norm_inv) * mu(i)
    end do
end subroutine laplacian






end module wawefunction
