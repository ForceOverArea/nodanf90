module solver1d
    use iso_fortran_env, only: dp => real64
    use generic1d
    use invert1d, only: invert_matrix_1d
    implicit none

contains

    function system_1d(ndof) result(sys)
        !! Initializes a system of `ndof` nodes to be linked with elements.
        !! The first element in the system will be initialized as the ground node.
        integer, intent(in) :: ndof
        type(node_1d), allocatable, dimension(:) :: sys

        allocate(sys(ndof))

        sys(1)%is_locked = .true. ! Ground the first node
        sys(1)%potential = 0
    end function system_1d

    function get_sys_dof(sys) result(dof)
        !! Gets an array of the degrees of freedom in a system
        type(node_1d), allocatable, dimension(:), intent(inout) :: sys
        type(node_1d), allocatable, dimension(:) :: dof
        integer :: i, ndof = 0

        do i = 1, size(sys)
            if (.not. sys(i)%is_locked) then
                ndof = ndof + 1
            end if
        end do

        allocate(dof(ndof))

        ndof = 1 ! reuse variable for the last index added to `dof`
        do i = 1, size(sys)
            if (.not. sys(i)%is_locked) then
                dof(i) = sys(i)
                ndof = ndof + 1
            end if
        end do
    end function get_sys_dof

    subroutine solve_1d(sys, margin, limit)
        !! Solves a nodal analysis problem
        type(node_1d), allocatable, dimension(:), intent(inout) :: sys
        real(dp), intent(in) :: margin
        integer, intent(in) :: limit
        type(node_1d), allocatable, dimension(:) :: dof
        real(dp), dimension(size(sys), 1) :: f_x, f_xdx
        real(dp), dimension(size(sys)) :: partial
        real(dp), dimension(size(sys), size(sys)) :: jacobian
        real(dp) :: dx = 1E-5
        integer :: ndof, i, count = 0

        dof = get_sys_dof(sys)
        ndof = size(dof)

        ! Calculate current error and prime loop
        f_x(:, 1) = flow_balances(sys)

        do while (sum(abs(f_x)) > margin)
            if (count > limit) then
                error stop 'despite reaching the iteration limit, no solution was found'
            end if
            do i = 1, ndof           
                ! Add `dx` to 1 potential at index `i`
                sys(i)%potential = sys(i)%potential + dx
                
                ! Calculate new error
                f_xdx(:, 1) = flow_balances(sys)
                
                ! Revert change by `dx`
                sys(i)%potential = sys(i)%potential - dx

                ! Calculate partial derivative and add to jacobian
                partial = (f_xdx - f_x(:, 1)) / dx
                jacobian(i, :) = partial(:) 
            end do

            ! f_xdx is reused here for storage for delta values
            f_xdx = invert_matrix_1d(jacobian) * f_x

            ! Adjust potentials by delta values
            do i = 1, ndof
                dof(i)%potential = dof(i)%potential - f_xdx(i, 1)
            end do

            ! Calculate current error for next iteration
            f_x(:, 1) = flow_balances(sys)

            ! Increment counter
            count = count + 1
        end do
    end subroutine solve_1d

end module solver1d