module test_functions_mod

    implicit none
    
contains

function rosenbrock(x) result(f)
    ! n-dimensional Rosenbrock function
    ! The global minimum is f(1,....,1) = 0

    implicit none

    real,dimension(:),intent(in) :: x
    real :: f

    integer :: N, i

    ! Get dimension
    N = size(x)

    ! Initialize
    f = 0.

    ! Loop through dimensions
    do i=1,N-1
        f = f + 100.*(x(i+1) - x(i)**2)**2 + (1. - x(i))**2
    end do
    
end function rosenbrock


function himmelblau(x) result(f)
    ! Himmelblau's function (2-dimensional)
    ! This function has 4 minima:
    ! f(3,2) = 0
    ! f(-2.805118,3.131312) = 0
    ! f(-3.779310,-3.283186) = 0
    ! f(3.584428,-1.848126) = 0

    implicit none
    
    real,dimension(:),intent(in) :: x
    real :: f

    f = (x(1)**2 + x(2) - 11.0)**2 + (x(1) + x(2)**2 - 7.)**2
    
end function himmelblau


end module test_functions_mod