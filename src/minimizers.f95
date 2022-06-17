module minimizers_mod

    implicit none

    type :: optimizer_result

        real,dimension(:),allocatable :: x_opt
        real :: f_opt
    
    end type optimizer_result
    
contains

function nelder_mead(f, x0) result(opt_result)
    ! Minimizes the function f

    implicit none
    
    interface
        real function f(x)
            real,dimension(:),allocatable,intent(in) :: x
        end function f
    end interface
    real,dimension(:),allocatable :: x0
    type(optimizer_result) :: opt_result

    integer :: N, i
    real,dimension(:,:),allocatable :: x_simp
    real,dimension(:),allocatable :: f_simp

    ! Get dimension
    N = size(x0)

    ! Initialize simplex vertices
    allocate(x_simp(N,N+1))
    x_simp(:,1) = x0
    do i=2,N+1
        x_simp(:,i) = x0
        x_simp(i-1,i) = x_simp(i-1,i) + 1.0
    end do

    ! Get function values for initial simplex
    allocate(f_simp(N+1))
    do i=1,N+1
        f_simp(i) = f(x_simp(:,i))
    end do

    
end function nelder_mead
    
end module minimizers_mod