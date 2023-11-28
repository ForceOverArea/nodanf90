module circuit1d
    use iso_fortran_env, only: dp => real64
    use generic1d, only: element_1d, node_1d
    implicit none

    ! 1D resistor element type
    type, extends(element_1d) :: resistor_1d
        real(dp) :: resistance = 1
        type(node_1d), pointer :: input, output
    contains
        procedure :: flow => resistor_1d_flow
    end type resistor_1d

    ! 1D battery element type
    type, extends(element_1d) :: battery_1d
        real(dp) :: voltage = 1
        type(node_1d), pointer :: input, output
        logical :: use_output = .false.
    contains
        procedure :: flow => battery_1d_flow
    end type battery_1d

contains

    ! Calculates current through a 1D resistor
    pure function resistor_1d_flow(this) result(current)
        class(resistor_1d), intent(in) :: this
        real(dp) :: current, voltage_drop

        if (this%resistance == 0) then
            error stop 'error: (procedure resistor_1d%flow): resistance cannot be zero!'
        end if

        voltage_drop = this%output%get_potential() - this%input%get_potential()
        current = voltage_drop / this%resistance
    end function resistor_1d_flow

    ! Calculates current through a 1D battery
    pure function battery_1d_flow(this) result(current)
        class(battery_1d), intent(in) :: this
        real(dp) :: current

        if (this%use_output) then
            current = this%output%flow_balance()
        else 
            current = this%input%flow_balance() 
        end if
    end function battery_1d_flow

end module circuit1d
