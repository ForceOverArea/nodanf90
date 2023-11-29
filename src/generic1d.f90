!=============================================
! Nodan Generics Module
! Grant Christiansen 2023
!=============================================
! Contains abstract 1D node and element types 
! for building more specific problem solvers
!=============================================
module generic1d
    use iso_fortran_env, only: dp => real64
    implicit none 

    ! Interface definition for a 1D element
    type, abstract :: element_1d
    contains
        procedure(flow_function), deferred :: flow
    end type element_1d

    abstract interface 
        pure function flow_function(this) result(absflow)
            import :: element_1d, dp
            class(element_1d), intent(in) :: this
            real(dp) :: absflow
        end function flow_function
    end interface

    ! Fixed-size type for allocatable arrays of elements
    type elem_box_1d
        class(element_1d), pointer :: elem
    end type elem_box_1d

    ! 1D Node type
    type node_1d
        real(dp) :: potential = 1, applied_step = 0
        logical :: is_locked
        type(elem_box_1d), allocatable, dimension(:) :: inputs, outputs
    contains
        procedure :: flow_balance    
        procedure :: get_potential
        procedure :: add_input
        procedure :: add_output
    end type node_1d

contains

    ! Elemental function for getting flows from nodal inputs/outputs
    elemental function get_flows(elems) result(flows)
        class(elem_box_1d), intent(in) :: elems
        real(dp) :: flows

        flows = elems%elem%flow()
    end function get_flows

    ! Procedure for performing nodal flow balance
    pure function flow_balance(this) result(discrepancy)
        class(node_1d), intent(in) :: this
        real(dp) :: discrepancy

        discrepancy = sum(get_flows(this%inputs)) - sum(get_flows(this%outputs))
    end function flow_balance

    elemental function flow_balances(nodes) result(discreps)
        type(node_1d), intent(in) :: nodes
        real(dp) :: discreps

        discreps = nodes%flow_balance()
    end function flow_balances

    ! Procedure for getting the potential of a node INCLUDING any
    ! step value applied by an attached element.
    pure function get_potential(this) result(potential)
        class(node_1d), intent(in) :: this
        real(dp) :: potential

        potential = this%potential + this%applied_step
    end function get_potential

    ! Adds an input to the nodal input list
    subroutine add_input(this, elem)
        class(node_1d), intent(inout) :: this
        class(element_1d), target, intent(in) :: elem
        type(elem_box_1d), dimension(size(this%inputs) + 1) :: tmp
        type(elem_box_1d), save :: newbox
        integer :: isize

        isize = size(this%inputs)

        newbox%elem => elem
        tmp(1:isize) = this%inputs(:)
        tmp(isize + 1) = newbox

        if (allocated(this%inputs)) then
            deallocate(this%inputs)
            allocate(this%inputs(isize + 1))
        else 
            allocate(this%inputs(1))
        end if

        this%inputs(:) = tmp(:)
    end subroutine add_input

    ! Adds an output to the nodal output list
    subroutine add_output(this, elem)
        class(node_1d), intent(inout) :: this
        class(element_1d), target, intent(in) :: elem
        type(elem_box_1d), dimension(size(this%outputs) + 1) :: tmp
        type(elem_box_1d), save :: newbox
        integer :: osize

        osize = size(this%outputs)

        newbox%elem => elem
        tmp(1:osize) = this%outputs(:)
        tmp(osize + 1) = newbox

        if (allocated(this%outputs)) then
            deallocate(this%outputs)
            allocate(this%outputs(osize + 1))
        else 
            allocate(this%outputs(1))
        end if

        this%outputs(:) = tmp(:)
    end subroutine add_output

end module generic1d