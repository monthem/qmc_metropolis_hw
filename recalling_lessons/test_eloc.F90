program test_eloc
    use H_local_energy
    implicit none

    double precision :: a, r(3)
    double precision :: E_loc

    r = [-2.0d0, 1.13d0, 5.0d0]
    a = 0.5d0
    E_loc = local_energy(a, r)
    print *, e_loc
end program



