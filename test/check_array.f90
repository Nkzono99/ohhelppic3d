module m_check_array
    use m_str
    use testdrive, only: error_type, test_failed
    implicit none

    private
    public check

    interface check
        module procedure :: check_int_array1d
        module procedure :: check_double_array3d
    end interface

contains

    subroutine check_int_array1d(error, arr, expected, message, more)
        type(error_type), allocatable, intent(out) :: error
        integer, intent(in) :: arr(:)
        integer, intent(in) :: expected(:)
        character(*), intent(in), optional :: message
        character(*), intent(in), optional :: more

        if (any(arr /= expected)) then
            if (present(message)) then
                call test_failed(error, message, more)
            else
                call test_failed(error, 'Integer array values missmatch', more)
            end if
        end if
    end subroutine

    subroutine check_double_array3d(error, arr, expected, th, message, more)
        type(error_type), allocatable, intent(out) :: error
        double precision, intent(in) :: arr(:, :, :)
        double precision, intent(in) :: expected(:, :, :)
        double precision, intent(in) :: th
        character(*), intent(in), optional :: message
        character(*), intent(in), optional :: more

        if (any(abs(arr - expected) > th)) then
            if (present(message)) then
                call test_failed(error, message, more)
            else
                call test_failed(error, 'Double array 3d values missmatch', more)
            end if
        end if
    end subroutine
end module
