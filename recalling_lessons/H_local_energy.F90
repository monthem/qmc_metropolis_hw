module H_local_energy
    implicit none
    contains

    function potential(r) result(V)
        ! V = -1 / r

        ! Argument
        double precision, intent(in) :: r(3)
        ! Output
        double precision :: V
        ! Local variable
        double precision :: norm_r

        norm_r = sqrt(sum(r**2))
        if (norm_r == 0.0) then
            print *, "Error: Division by zero in potential calculation!"
            stop
        else
            V = -1.0d0 / norm_r
        end if
    end function potential

    function kinetic(a, r) result(T)
        ! Calculates the local kinetic energy T_l(r) = -0.5 (a^2 - 2a/r); derived by applying the Laplacian operator to the wavefunction
        
        ! Arguments
        double precision, intent(in) :: a, r(3)
        ! Output
        double precision             :: T
        ! Local variables
        double precision             :: norm_r

        norm_r = sqrt(sum(r**2))
        if (norm_r == 0) stop "Error: Division by zero in kinetic energy calculation!"

        T = -0.5d0 * (a**2 - 2.0d0 * a / norm_r)

    end function kinetic

    function wawefunction(a, r) result(psi)
        ! A function that evaluates psi(r) at a given r, psi = exp(-ar)
        
        ! Arguments
        double precision, intent(in) :: a, r(3)
        ! Output
        double precision :: psi

        psi = exp(-a*sqrt(sum(r**2)))
    end function wawefunction 

    function local_energy(a, r) result(E_loc)
        ! Returns E_loc = T + V
        
        ! Arguments
        double precision, intent(in)        :: a, r(3)
        ! Output
        double precision                    :: E_loc

        E_loc = kinetic(a, r) + potential(r)
    end function local_energy
    
end module H_local_energy
