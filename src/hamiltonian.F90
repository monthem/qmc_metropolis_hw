module hamiltonian
use wawefunction
implicit none
contains

function nuc_nuc(A, Z, R) result(V_NN)
    ! Arguments 
    integer, intent(in)             :: A, Z(A)
    double precision, intent(in)    :: R(A, 3)
    ! Output
    double precision                :: V_NN
    ! Local variable
    integer                         :: i, j
    double precision                :: norm_R_AB

    V_NN = 0d0
    do i = 1, A
        do j = i + 1, A
            print *, "computing atom ", i, "with atom", j
            norm_R_AB = sqrt(sum((R(i,:) - R(j,:))**2))
            print *, "norm is: ", norm_R_AB
            if (norm_R_AB <= 0.5d0) stop "Attempted nuclear fusion in function nuc_nuc!" ! Don't get two H atoms closer than 0.5 Bohr to each other
            V_NN = V_NN + Z(i) * Z(j) / norm_R_AB
            print *, Z(i) * Z(j) / norm_R_AB
        end do
    end do
end function nuc_nuc

function nuc_el(A, Z, q) result (V_NE)
    ! NOT A TOTAL POTENTIAL BUT JUST THAT OF ONE ELECTRON
    ! Arguments
    integer, intent(in)             :: A, Z(A)
    double precision, intent(in)    :: q(A, 3)
    ! Output
    double precision                :: V_NE
    ! Local variable 
    integer                         :: i
    double precision                :: norm_q

    V_NE = 0d0
    do i = 1, A
        norm_q = sqrt(sum(q(i,:)**2))
        if (norm_q <= 1d-4) stop "Electron too close to a nucleus in function nuc_el!"
        V_NE = V_NE - Z(i)/norm_q
    end do
end function nuc_el

function el_el(n, A, big_Q) result(V_EE)
    ! Arguments
    integer                 :: n, A
    double precision        :: big_Q(n, A, 3)
    double precision        :: V_EE
    ! Local variable
    integer                 :: i, j
    double precision        :: r_ij(3), norm_r
    
     V_EE = 0d0
    do i = 1, n
        do j = i + 1, n
            r_ij = big_Q(i, 1, :) - big_Q(j, 1, :) ! Positions of electrons both relative to nuclei 1
            norm_r = sqrt(sum(r_ij**2))
            if (norm_r <= 1d-6) stop "Electrons too close to each other in function el_el!)"
            V_EE = V_EE + 1d0 / norm_r
        end do
    end do
end function el_el

function kinetic_energy(A, coeffs, alpha, q, mu, phi) result(T)
    ! Arguments
    integer, intent(in)             :: A
    double precision, intent(in)    :: coeffs(A), alpha, q(A, 3), mu(A), phi
    ! Output
    double precision                :: T
    ! Local variables
    double precision                :: del_sq
    
    del_sq = 0d0
    call laplacian(A, coeffs, alpha, q, mu, del_sq)
    T = -0.5d0 * del_sq / phi
end function kinetic_energy



end module hamiltonian
        



