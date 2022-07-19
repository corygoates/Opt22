module minimizers_mod

    use sort_mod

    implicit none

    type :: optimizer_result

        real,dimension(:),allocatable :: x_opt
        real :: f_opt
    
    end type optimizer_result
    
contains

function nelder_mead(f, x0, termination_tol, reflection_coef, expansion_coef, contraction_coef, shrink_coef) result(opt_result)
    ! Minimizes the function f

    implicit none
    
    interface
        real function f(x)
            real,dimension(:),intent(in) :: x
        end function f
    end interface
    real,dimension(:),allocatable :: x0
    real,optional :: termination_tol, reflection_coef, expansion_coef, contraction_coef, shrink_coef
    type(optimizer_result) :: opt_result

    integer :: N, i, iteration
    real :: avg, std_dev
    real :: tol, alpha, gamma, rho, sigma, f_r, f_e, f_c
    real,dimension(:,:),allocatable :: x_simp
    real,dimension(:),allocatable :: f_simp, x_cent, x_r, x_e, x_c
    integer,dimension(:),allocatable :: i_sorted

    ! Set defaults
    if (present(termination_tol)) then
        tol = termination_tol
    else
        tol = 1.e-6
    end if

    if (present(reflection_coef)) then
        alpha = reflection_coef
    else
        alpha = 1.
    end if

    if (present(expansion_coef)) then
        gamma = expansion_coef
    else
        gamma = 2.
    end if

    if (present(contraction_coef)) then
        rho = contraction_coef
    else
        rho = 0.5
    end if

    if (present(shrink_coef)) then
        sigma = shrink_coef
    else
        sigma = 0.5
    end if

    ! Get dimension
    N = size(x0)

    ! Initialize simplex vertices
    allocate(x_simp(N,N+1))
    allocate(x_cent(N))
    allocate(x_r(N))
    allocate(x_e(N))
    allocate(x_c(N))
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

    ! Enter iteration
    iteration = 0
    do while (.true.)

        ! Update iteration number
        iteration = iteration + 1
        write(*,*)
        write(*,*)
        write(*,*) "Iteration: ", iteration
        write(*,*) "x:"
        do i=1,N+1
            write(*,*) x_simp(:,i)
        end do
        write(*,*)
        write(*,*) "f: ", f_simp

        ! Check standard deviation for termination
        avg = sum(f_simp)/(N+1)
        std_dev = sum((f_simp-avg)**2)/N
        if (std_dev <= 1.e-12) exit

        ! Arrange the vertices
        call insertion_arg_sort(f_simp, i_sorted)

        ! Calculate centroid
        x_cent = x_simp(:,1)
        do i=2,N+1
            x_cent = x_cent + x_simp(:,i)
        end do
        x_cent = x_cent / (N+1)

        ! Reflection
        write(*,*)
        write(*,*) "Reflecting"
        x_r = x_cent + alpha*(x_cent - x_simp(:,i_sorted(N+1)))
        f_r = f(x_r)
        write(*,*) "x_r: ", x_r
        write(*,*) "f_r: ", f_r

        ! Check if the reflected point is good but not too good
        if (f_r >= f_simp(i_sorted(1)) .and. f_r < f_simp(i_sorted(N))) then
            x_simp(:,i_sorted(N+1)) = x_r
            f_simp(i_sorted(N+1)) = f_r
            cycle

        ! Expansion
        else if (f_r < f_simp(i_sorted(1))) then
            write(*,*)
            write(*,*) "Expanding"

            ! Calculate expanded point
            x_e = x_cent + gamma*(x_r - x_cent)
            f_e = f(x_e)
            write(*,*) "x_e: ", x_e
            write(*,*) "f_e: ", f_e

            ! Compare to reflected point
            if (f_e < f_r) then
                x_simp(:,i_sorted(N+1)) = x_e
                f_simp(i_sorted(N+1)) = f_e
            else
                x_simp(:,i_sorted(N+1)) = x_r
                f_simp(i_sorted(N+1)) = f_r
            end if

            cycle

        ! Contraction
        else

            write(*,*)
            write(*,*) "Contracting"

            ! Check if we're better than the worst
            if (f_r < f_simp(i_sorted(N+1))) then

                ! Calculate contracted point
                x_c = x_cent + rho*(x_r - x_cent)
                f_c = f(x_c)
                write(*,*) "x_c: ", x_c
                write(*,*) "f_c: ", f_c

                ! Compare to reflected point
                if (f_c < f_r) then
                    x_simp(:,i_sorted(N+1)) = x_c
                    f_simp(i_sorted(N+1)) = f_c
                    cycle
                end if

            ! We're worse than the worst
            else

                ! Calculate contracted point
                x_c = x_cent + rho*(x_simp(:,i_sorted(N+1)) - x_cent)
                f_c = f(x_c)
                write(*,*) "x_c: ", x_c
                write(*,*) "f_c: ", f_c

                ! Compare to worst point
                if (f_c < f_simp(i_sorted(N+1))) then
                    x_simp(:,i_sorted(N+1)) = x_c
                    f_simp(i_sorted(N+1)) = f_c
                    cycle
                end if

            end if

        end if
        
        ! Shrink
        write(*,*)
        write(*,*) "Shrinking"
        ! If we've made it here, then the simplex needs to shrink
        do i=2,N+1

            ! Shrink in vertices
            x_simp(:,i_sorted(i)) = x_simp(:,i_sorted(1)) + sigma*(x_simp(:,i_sorted(i)) - x_simp(:,i_sorted(1)))

            ! Calculate new values
            f_simp(i_sorted(i)) = f(x_simp(:,i_sorted(i)))

        end do

    end do

    
end function nelder_mead
    
end module minimizers_mod