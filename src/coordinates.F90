module coordinates
implicit none
contains

subroutine nuclei2bohr(A, nuclei)
    ! Arguments
    integer, intent(in)                         :: A
    double precision,intent(inout)              :: nuclei(A, 3)
    ! Local argument
    double precision, parameter                 :: ang2bohr = 1.8897259886d0
    
    nuclei = nuclei * ang2bohr
end subroutine nuclei2bohr

subroutine form_Q(A, n, nuclei, electrons, Q)
    ! Arguments
    integer, intent(in)             :: A, n
    double precision, intent(in)    :: nuclei(A, 3), electrons(n, 3)
    double precision, intent(out)   :: Q(n, A, 3)
    ! Local variables
    integer                         :: i, j, k

    do i = 1, n                     ! Loop over each electron coordinate
        do j = 1, A                 ! Loop over each nuclear coordinate
            Q(i, j, :) = electrons(i, :) - nuclei(j, :)
        end do
    end do
end subroutine form_Q

end module coordinates
