module wawefunction     ! functions for the calculation of individual AOs mu, MOs phi and the total WF psi
implicit none
contains

subroutine AO(A, coeffs, alphas, q, mu)
    ! subroutine that takes the number of nuclei A, an array of AO coefficients and exponents of size N, an array of position vectors q of size N by 3, and returns an array of values of the AOs mu of size N
    ! Arguments
    integer, intent(in)             :: A 
    double precision, intent(in)    :: coeffs(A), alphas(A), q(A, 3)
    double precision, intent(out)   :: mu(A)
    ! Local variable
    integer                         :: i

    do i = 1, A
        mu(i) = coeffs(i) * exp(-alphas(i)*sqrt(sum(q(i,:)**2)))
    end do
end subroutine AO

subroutine MO(A, coeffs, alphas, q, mu, phi)
    ! Arguments
    integer, intent(in)                 :: A
    double precision, intent(in)        :: coeffs(A), alphas(A), q(A, 3)
    double precision, intent(inout)     :: mu(A)
    double precision, intent(out)       :: phi
    ! Local variable
    integer                             :: i

    call AO(A, coeffs, alphas, q, mu)
    phi = sum(mu)
end subroutine MO

function WF(MO_1, MO_2) result(psi)
    ! Arguments
    double precision, intent(in)        :: MO_1, MO_2
    ! Output 
    double precision                    :: psi
    psi = MO_1 * MO_2
end function WF

end module wawefunction
