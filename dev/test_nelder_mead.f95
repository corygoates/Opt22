program test_nelder_mead

    use test_functions_mod
    use minimizers_mod

    implicit none

    real,dimension(:),allocatable :: x0
    type(optimizer_result) :: solution

    ! Set up initial guess
    allocate(x0(2))
    x0(1) = -1.
    x0(2) = 3.

    ! Run optimizer
    solution = nelder_mead(rosenbrock, x0)
    
end program test_nelder_mead