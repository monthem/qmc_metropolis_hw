program test_wf
    use H_local_energy
    implicit none

    double precision :: T, a, r(3)

    a = 0.5d0
    r = [-1.0d0, 0.88d0, -0.542d0]

    T = kinetic(a, r)

    print *, T

    r = [0.0d0, 0.0d0, 0.0d0]
    T = kinetic(a, r)
    print *, T
end program test_wf
