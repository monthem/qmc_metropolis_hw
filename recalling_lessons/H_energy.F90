module H_energy
    use H_local_energy
    implicit none
    contains

    subroutine numerical_estimate
        ! Performs a numerical integration of a H atom wf on a below specified grid
        ! E_el = integral of |psi(r)|²E_loc dr over integral of |psi(r)|²dr
        
        double precision            :: dx = 10.0d0/49.0d0, dy, dz ! (5+5)/(50-1)
        integer                     :: i, j, k
        double precision            :: a = 1.05d0, x, y, z

        double precision            :: E_mean, E_i, w_i, r(3)
        double precision            :: numerator = 0, denominator = 0
        
        ! Variance calculation variables
        double precision            :: sigma_sq
        double precision            :: sq_mean = 0.0d0 , mean_sq = 0.0d0 ! sq_mean is <E_loc²>, mean_sq is <E_loc>^2

        dy = dx
        dz = dx

        ! Loop generates values of x, y, z from -5 to 5 in steps of ~0.2041 (avoids having x, y, z = 0)
        do i = 1, 50 
            x = -5.0d0 + (i - 1) * dx
            do j = 1, 50
                y = -5.0d0 + (j - 1) * dx
                do k = 1, 50
                    z = -5.0d0 + (k - 1) * dx

                    r = [x, y, z]
 
                    E_i = local_energy(a, r)
                    sq_mean = sq_mean + E_i**2 ! sums the E_i^2 to perform the mean later
                    mean_sq = mean_sq + E_i ! sums E_i to perform the mean later and square it
                    w_i = wawefunction(a, r)**2 * dx*dy*dz

                    numerator = numerator + w_i * E_i
                    denominator = denominator + w_i
                end do
            end do
        end do
        E_mean = numerator / denominator
        sq_mean = sq_mean / 50.0d0**3 ! 50**3 is because there is 50**3 evaluations of E_i
        mean_sq = (mean_sq / 50.0d0**3)**2 
        print *, "E = ", E_mean, " a.u"
        print *, "variance of the local energy is: ", sq_mean - mean_sq
    end subroutine numerical_estimate

    
end module H_energy
