program test_nelder_mead

    use test_functions_mod
    use minimizers_mod

    implicit none

    real,dimension(:),allocatable :: x0
    type(optimizer_result) :: solution

    ! Set up initial guess
    allocate(x0(2))
    x0(1) = -100.
    x0(2) = 300.

    ! Run optimizer
    solution = nelder_mead(rosenbrock, x0, termination_tol=1e-14)
    write(*,*)
    write(*,*) "x_opt: ", solution%x_opt
    write(*,*) "f_opt: ", solution%f_opt
    write(*,*) "Iterations: ", solution%iterations
    
end program test_nelder_mead