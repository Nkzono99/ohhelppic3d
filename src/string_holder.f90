module m_string_holder
    implicit none

    private
    public t_StringHolder
    public new_StringHolder

    type :: t_StringHolder
        character(len=:), allocatable :: string
    end type

contains

    function new_StringHolder(string) result(obj)
        character(*), intent(in) :: string
        type(t_StringHolder) :: obj

        obj%string = string
    end function

end module
