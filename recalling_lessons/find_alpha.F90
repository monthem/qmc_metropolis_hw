program find_optimal_alpha
    use H_local_energy
    implicit none
    
    integer, parameter          :: n = 30
    double precision            :: alphas(6), vectors(n, 3)
    double precision            :: a, r(3), x, dx, E_l
    integer                     :: i, j

    ! generates 30 vectors with x = -5, ..., +5 and != 0
    dx = 10.0d0 / (n - 1)
    do i = 1, n
        x = -5.0d0 + (i - 1) * dx ! go from -5 to +5 in 30 steps
        vectors(i, :) = [x, 0.0d0, 0.0d0]
    end do
    alphas = [0.1d0, 0.2d0, 0.5d0, 1.0d0, 1.5d0, 2.0d0]

    do i = 1, n
        x = vectors(i, 1)
        r = vectors(i, :)

        write(10, '(F5.2)', advance="no") x

        do j = 1, 6
            a = alphas(j)
            E_l = local_energy(a, r)
            write(10, '(2X, F10.4)', advance="no") E_l
        end do
        
        write(10, *)
    end do

end program find_optimal_alpha

