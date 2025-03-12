program test_potential
    use H_local_energy
    implicit none

    double precision :: r(3)
    double precision :: V

    r = [-2.0d0, 1.13d0, 5.0d0]

    V = potential(r)

    print *, V

    r = [0.0d0, 0.0d0, 0.0d0]
    V = potential(r)
    print *, V
end program



