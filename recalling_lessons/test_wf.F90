program test_wf
    use H_local_energy
    implicit none

    double precision :: a, r(3), psi

    a = 0.5d0
    r = [-1.0d0, 0.88d0, -0.542d0]

    psi = wawefunction(a, r)

    print *, psi
end program test_wf
