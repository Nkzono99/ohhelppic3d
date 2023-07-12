module m_toml_wrapper
    use mpi
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    use tomlf, only: toml_table, &
                     toml_array, &
                     toml_parse, &
                     get_value, &
                     toml_serialize, &
                     len, toml_error, toml_stat
    use m_string_holder

    implicit none

    private
    public t_TomlWrapper

    type :: t_TomlWrapper
        type(toml_table), allocatable :: table

    contains

        procedure :: load => toml_load
        procedure :: to_string => toml_to_string

        generic :: require_table => &
            require_table1, require_table2, require_table3
        procedure, private :: require_table1 => toml_require_table1
        procedure, private :: require_table2 => toml_require_table2
        procedure, private :: require_table3 => toml_require_table3

        generic :: require_array => &
            require_array1, require_array2, require_array3, require_array4
        procedure, private :: require_array1 => toml_require_array1
        procedure, private :: require_array2 => toml_require_array2
        procedure, private :: require_array3 => toml_require_array3
        procedure, private :: require_array4 => toml_require_array4

        generic :: require_int => require_int1, require_int2, require_int3, require_int4
        procedure, private :: require_int1 => toml_require_int1
        procedure, private :: require_int2 => toml_require_int2
        procedure, private :: require_int3 => toml_require_int3
        procedure, private :: require_int4 => toml_require_int4

        generic :: require_real => &
            require_real1, require_real2, require_real3, require_real4
        procedure, private :: require_real1 => toml_require_real1
        procedure, private :: require_real2 => toml_require_real2
        procedure, private :: require_real3 => toml_require_real3
        procedure, private :: require_real4 => toml_require_real4

        generic :: require_double => &
            require_double1, require_double2, require_double3, require_double4
        procedure, private :: require_double1 => toml_require_double1
        procedure, private :: require_double2 => toml_require_double2
        procedure, private :: require_double3 => toml_require_double3
        procedure, private :: require_double4 => toml_require_double4

        procedure :: require_string1 => toml_require_string1
        procedure :: require_string2 => toml_require_string2
        procedure :: require_string3 => toml_require_string3
        procedure :: require_string4 => toml_require_string4

        generic :: require_int_array => &
            require_int_array1, require_int_array2, require_int_array3, require_int_array4
        procedure, private :: require_int_array1 => toml_require_int_array1
        procedure, private :: require_int_array2 => toml_require_int_array2
        procedure, private :: require_int_array3 => toml_require_int_array3
        procedure, private :: require_int_array4 => toml_require_int_array4

        generic :: require_real_array => &
            require_real_array1, require_real_array2, require_real_array3, require_real_array4
        procedure, private :: require_real_array1 => toml_require_real_array1
        procedure, private :: require_real_array2 => toml_require_real_array2
        procedure, private :: require_real_array3 => toml_require_real_array3
        procedure, private :: require_real_array4 => toml_require_real_array4

        generic :: require_double_array => &
            require_double_array1, require_double_array2, require_double_array3, require_double_array4
        procedure, private :: require_double_array1 => toml_require_double_array1
        procedure, private :: require_double_array2 => toml_require_double_array2
        procedure, private :: require_double_array3 => toml_require_double_array3
        procedure, private :: require_double_array4 => toml_require_double_array4

        generic :: require_string_array => &
            require_string_array1, require_string_array2, require_string_array3, require_string_array4
        procedure, private :: require_string_array1 => toml_require_string_array1
        procedure, private :: require_string_array2 => toml_require_string_array2
        procedure, private :: require_string_array3 => toml_require_string_array3
        procedure, private :: require_string_array4 => toml_require_string_array4

        generic :: require_string_array2d => &
            require_string_array2d1, require_string_array2d2, require_string_array2d3, require_string_array2d4
        procedure, private :: require_string_array2d1 => toml_require_string_array2d1
        procedure, private :: require_string_array2d2 => toml_require_string_array2d2
        procedure, private :: require_string_array2d3 => toml_require_string_array2d3
        procedure, private :: require_string_array2d4 => toml_require_string_array2d4
    end type

contains

    subroutine error_stop(message)
        character(len=*), intent(in) :: message
        integer :: ierr

        write (stderr, '(a)') 'Error: '//message
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        stop 1
    end subroutine

    subroutine toml_load(self, toml_filepath)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: toml_filepath

        integer :: io
        type(toml_error), allocatable :: error

        open (file=toml_filepath, newunit=io, status='old')
        call toml_parse(self%table, io, error)
        close (io)

        if (allocated(error)) then
            call error_stop(error%message)
        end if
    end subroutine

    function toml_require_table1(self, name) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name
        type(toml_table), pointer :: ret

        integer :: stat

        call get_value(self%table, name, ret, stat=stat, requested=.false.)
        if (.not. associated(ret)) then
            call error_stop('Table '//name//' is not found')
        end if
    end function

    function toml_require_table2(self, name1, name2) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        type(toml_table), pointer :: ret

        integer :: stat

        type(toml_table), pointer :: table

        table => self%require_table(name1)

        call get_value(table, name2, ret, stat=stat, requested=.false.)

        if (.not. associated(ret)) then
            call error_stop('Table '//name1//'.'//name2//' is not found')
        end if
    end function

    function toml_require_table3(self, name1, name2, name3) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        type(toml_table), pointer :: ret

        integer :: stat

        call get_value(self%require_table(name1, name2), name3, ret, &
                       stat=stat, requested=.false.)

        if (.not. associated(ret)) then
            call error_stop('Table '//name1//'.'//name2//'.'//name3//' is not found')
        end if
    end function

    function toml_require_array1(self, name) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name
        type(toml_array), pointer :: ret

        integer :: stat

        call get_value(self%table, name, ret, stat=stat)
        if (stat /= toml_stat%success) then
            call error_stop('Array '//name//' is not found')
        end if
    end function

    function toml_require_array2(self, name1, name2) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        type(toml_array), pointer :: ret

        integer :: stat

        call get_value(self%require_table(name1), name2, ret, stat=stat)
        if (stat /= toml_stat%success) then
            call error_stop('Array '//name1//','//name2//' is not found')
        end if
    end function

    function toml_require_array3(self, name1, name2, name3) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        type(toml_array), pointer :: ret

        integer :: stat

        call get_value(self%require_table(name1, name2), name3, ret, stat=stat)
        if (stat /= toml_stat%success) then
            call error_stop('Array '//name1//','//name2//','//name3//' is not found')
        end if
    end function

    function toml_require_array4(self, name1, name2, name3, name4) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        character(len=*), intent(in) :: name4
        type(toml_array), pointer :: ret

        integer :: stat

        call get_value(self%require_table(name1, name2, name3), name4, ret, stat=stat)
        if (stat /= toml_stat%success) then
            call error_stop('Array '//name1//','//name2//','//name3//','//name4//' is not found')
        end if
    end function

    function toml_require_int1(self, name, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer, optional, intent(in) :: default
        integer :: ret

        integer :: stat

        call get_value(self%table, name, ret, default, stat=stat)
        if (stat /= toml_stat%success) then
            call error_stop('Integer '//name//' is not found')
        end if
    end function

    function toml_require_int2(self, name1, name2, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        integer, optional, intent(in) :: default
        integer :: ret

        integer :: stat

        call get_value(self%require_table(name1), name2, ret, default, stat=stat)

        if (stat /= toml_stat%success) then
            call error_stop('Integer '//name1//'.'//name2//' is not found')
        end if
    end function

    function toml_require_int3(self, name1, name2, name3, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        integer, optional, intent(in) :: default
        integer :: ret

        integer :: stat

        call get_value(self%require_table(name1, name2), name3, ret, default, stat=stat)

        if (stat /= toml_stat%success) then
            call error_stop('Integer '//name1//'.'//name2//','//name3//' is not found')
        end if
    end function

    function toml_require_int4(self, name1, name2, name3, name4, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        character(len=*), intent(in) :: name4
        integer, optional, intent(in) :: default
        integer :: ret

        integer :: stat

        call get_value(self%require_table(name1, name2, name3), name4, ret, default, stat=stat)

        if (stat /= toml_stat%success) then
            call error_stop('Integer '//name1//'.'//name2//','//name3//','//name4//' is not found')
        end if
    end function

    function toml_require_real1(self, name, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name
        real, optional, intent(in) :: default
        real  :: ret

        integer :: stat

        call get_value(self%table, name, ret, default, stat=stat)
        if (stat /= toml_stat%success) then
            call error_stop('Real '//name//' is not found')
        end if
    end function

    function toml_require_real2(self, name1, name2, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        real, optional, intent(in) :: default
        real  :: ret

        integer :: stat

        call get_value(self%require_table(name1), name2, ret, default, stat=stat)

        if (stat /= toml_stat%success) then
            call error_stop('Real '//name1//'.'//name2//' is not found')
        end if
    end function

    function toml_require_real3(self, name1, name2, name3, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        real, optional, intent(in) :: default
        real  :: ret

        integer :: stat

        call get_value(self%require_table(name1, name2), name3, ret, default, stat=stat)

        if (stat /= toml_stat%success) then
            call error_stop('Real '//name1//'.'//name2//','//name3//' is not found')
        end if
    end function

    function toml_require_real4(self, name1, name2, name3, name4, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        character(len=*), intent(in) :: name4
        real, optional, intent(in) :: default
        real :: ret

        integer :: stat

        call get_value(self%require_table(name1, name2, name3), name4, ret, default, stat=stat)

        if (stat /= toml_stat%success) then
            call error_stop('Real '//name1//'.'//name2//','//name3//','//name4//' is not found')
        end if
    end function

    function toml_require_double1(self, name, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name
        double precision, optional, intent(in) :: default
        double precision :: ret

        integer :: stat

        call get_value(self%table, name, ret, default, stat=stat)
        if (stat /= toml_stat%success) then
            call error_stop('Double precision '//name//' is not found')
        end if
    end function

    function toml_require_double2(self, name1, name2, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        double precision, optional, intent(in) :: default
        double precision :: ret

        integer :: stat

        call get_value(self%require_table(name1), name2, ret, default, stat=stat)

        if (stat /= toml_stat%success) then
            call error_stop('Double precision '//name1//'.'//name2//' is not found')
        end if
    end function

    function toml_require_double3(self, name1, name2, name3, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        double precision, optional, intent(in) :: default
        double precision :: ret

        integer :: stat

        call get_value(self%require_table(name1, name2), name3, ret, default, stat=stat)

        if (stat /= toml_stat%success) then
            call error_stop('Double precision '//name1//'.'//name2//','//name3//' is not found')
        end if
    end function

    function toml_require_double4(self, name1, name2, name3, name4, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        character(len=*), intent(in) :: name4
        double precision, optional, intent(in) :: default
        double precision :: ret

        integer :: stat

        call get_value(self%require_table(name1, name2, name3), name4, ret, default, stat=stat)

        if (stat /= toml_stat%success) then
            call error_stop('Double precision '//name1//'.'//name2//','//name3//','//name4//' is not found')
        end if
    end function

    function toml_require_string1(self, name, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name
        character(len=*), optional, intent(in) :: default
        character(len=:), allocatable :: ret

        integer :: stat

        call get_value(self%table, name, ret, default, stat=stat)
        if (stat /= toml_stat%success) then
            call error_stop('String '//name//' is not found')
        end if
    end function

    function toml_require_string2(self, name1, name2, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), optional, intent(in) :: default
        character(len=:), allocatable :: ret

        integer :: stat

        call get_value(self%require_table(name1), name2, ret, default, stat=stat)

        if (stat /= toml_stat%success) then
            call error_stop('String '//name1//'.'//name2//' is not found')
        end if
    end function

    function toml_require_string3(self, name1, name2, name3, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        character(len=*), optional, intent(in) :: default
        character(len=:), allocatable :: ret

        integer :: stat

        call get_value(self%require_table(name1, name2), name3, ret, default, stat=stat)

        if (stat /= toml_stat%success) then
            call error_stop('String '//name1//'.'//name2//','//name3//' is not found')
        end if
    end function

    function toml_require_string4(self, name1, name2, name3, name4, default) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        character(len=*), intent(in) :: name4
        character(len=:), allocatable, optional, intent(in) :: default
        character(len=:), allocatable :: ret

        integer :: stat

        call get_value(self%require_table(name1, name2, name3), name4, ret, default, stat=stat)

        if (stat /= toml_stat%success) then
            call error_stop('String '//name1//'.'//name2//','//name3//','//name4//' is not found')
        end if
    end function

    function toml_require_int_array1(self, name) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer, allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival), stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('Integer '//name//' is not found')
            end if
        end do
    end function

    function toml_require_int_array2(self, name1, name2) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        integer, allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name1, name2)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival), stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('Integer '//name1//'.'//name2//' is not found')
            end if
        end do
    end function

    function toml_require_int_array3(self, name1, name2, name3) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        integer, allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name1, name2, name3)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival), stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('Integer '//name1//'.'//name2//'.'//name3//' is not found')
            end if
        end do
    end function

    function toml_require_int_array4(self, name1, name2, name3, name4) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        character(len=*), intent(in) :: name4
        integer, allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name1, name2, name3, name4)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival), stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('Integer '//name1//'.'//name2//'.'//name3//'.'//name4//' is not found')
            end if
        end do
    end function

    function toml_require_real_array1(self, name) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name
        real, allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival), stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('Real '//name//' is not found')
            end if
        end do
    end function

    function toml_require_real_array2(self, name1, name2) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        real, allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name1, name2)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival), stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('Real '//name1//'.'//name2//' is not found')
            end if
        end do
    end function

    function toml_require_real_array3(self, name1, name2, name3) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        real, allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name1, name2, name3)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival), stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('Real '//name1//'.'//name2//'.'//name3//' is not found')
            end if
        end do
    end function

    function toml_require_real_array4(self, name1, name2, name3, name4) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        character(len=*), intent(in) :: name4
        real, allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name1, name2, name3, name4)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival), stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('Real '//name1//'.'//name2//'.'//name3//'.'//name4//' is not found')
            end if
        end do
    end function

    function toml_require_double_array1(self, name) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name
        double precision, allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival), stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('Double precision '//name//' is not found')
            end if
        end do
    end function

    function toml_require_double_array2(self, name1, name2) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        double precision, allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name1, name2)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival), stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('Double precision '//name1//'.'//name2//' is not found')
            end if
        end do
    end function

    function toml_require_double_array3(self, name1, name2, name3) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        double precision, allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name1, name2, name3)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival), stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('Double precision '//name1//'.'//name2//'.'//name3//' is not found')
            end if
        end do
    end function

    function toml_require_double_array4(self, name1, name2, name3, name4) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        character(len=*), intent(in) :: name4
        double precision, allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name1, name2, name3, name4)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival), stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('Double precision '//name1//'.'//name2//'.'//name3//'.'//name4//' is not found')
            end if
        end do
    end function

    function toml_require_string_array1(self, name) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name
        type(t_StringHolder), allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival)%string, stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('String array '//name//' is not found')
            end if
        end do
    end function

    function toml_require_string_array2(self, name1, name2) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        type(t_StringHolder), allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name1, name2)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival)%string, stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('String array '//name1//'.'//name2//' is not found')
            end if
        end do
    end function

    function toml_require_string_array3(self, name1, name2, name3) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        type(t_StringHolder), allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name1, name2, name3)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival)%string, stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('String array '//name1//'.'//name2//'.'//name3//' is not found')
            end if
        end do
    end function

    function toml_require_string_array4(self, name1, name2, name3, name4) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        character(len=*), intent(in) :: name4
        type(t_StringHolder), allocatable :: ret(:)

        integer :: stat

        type(toml_array), pointer :: array
        integer :: ival

        array => self%require_array(name1, name2, name3, name4)

        allocate (ret(len(array)))
        do ival = 1, size(ret)
            call get_value(array, ival, ret(ival)%string, stat=stat)

            if (stat /= toml_stat%success) then
                call error_stop('String array '//name1//'.'//name2//'.'//name3//'.'//name4//' is not found')
            end if
        end do
    end function

    function toml_require_string_array2d1(self, name) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name
        type(t_StringHolder), allocatable :: ret(:, :)

        integer :: stat

        type(toml_array), pointer :: array1
        type(toml_array), pointer :: array2
        integer :: ival1, ival2

        array1 => self%require_array(name)

        do ival1 = 1, len(array1)
            call get_value(array1, ival1, array2, stat=stat)
            if (.not. allocated(ret)) then
                allocate (ret(len(array2), len(array1)))
            end if

            if (stat /= toml_stat%success) then
                call error_stop('String array 2d '//name//' is not found')
            end if

            do ival2 = 1, len(array2)
                call get_value(array2, ival1, ret(ival1, ival2)%string, stat=stat)

                if (stat /= toml_stat%success) then
                    call error_stop('String array 2d '//name//' is not found')
                end if
            end do
        end do
    end function

    function toml_require_string_array2d2(self, name1, name2) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        type(t_StringHolder), allocatable :: ret(:, :)

        integer :: stat

        type(toml_array), pointer :: array1
        type(toml_array), pointer :: array2
        integer :: ival1, ival2

        array1 => self%require_array(name1, name2)

        do ival1 = 1, len(array1)
            call get_value(array1, ival1, array2, stat=stat)
            if (.not. allocated(ret)) then
                allocate (ret(len(array2), len(array1)))
            end if

            if (stat /= toml_stat%success) then
                call error_stop('String array 2d '//name1//'.'//name2//' is not found')
            end if

            do ival2 = 1, len(array2)
                call get_value(array2, ival1, ret(ival1, ival2)%string, stat=stat)

                if (stat /= toml_stat%success) then
                    call error_stop('String array 2d '//name1//'.'//name2//' is not found')
                end if
            end do
        end do
    end function

    function toml_require_string_array2d3(self, name1, name2, name3) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        type(t_StringHolder), allocatable :: ret(:, :)

        integer :: stat

        type(toml_array), pointer :: array1
        type(toml_array), pointer :: array2
        integer :: ival1, ival2

        array1 => self%require_array(name1, name2, name3)

        do ival1 = 1, len(array1)
            call get_value(array1, ival1, array2, stat=stat)
            if (.not. allocated(ret)) then
                allocate (ret(len(array2), len(array1)))
            end if

            if (stat /= toml_stat%success) then
                call error_stop('String array 2d '//name1//'.'//name2//' is not found')
            end if

            do ival2 = 1, len(array2)
                call get_value(array2, ival1, ret(ival1, ival2)%string, stat=stat)

                if (stat /= toml_stat%success) then
                    call error_stop('String array 2d '//name1//'.'//name2//' is not found')
                end if
            end do
        end do
    end function

    function toml_require_string_array2d4(self, name1, name2, name3, name4) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=*), intent(in) :: name1
        character(len=*), intent(in) :: name2
        character(len=*), intent(in) :: name3
        character(len=*), intent(in) :: name4
        type(t_StringHolder), allocatable :: ret(:, :)

        integer :: stat

        type(toml_array), pointer :: array1
        type(toml_array), pointer :: array2
        integer :: ival1, ival2

        array1 => self%require_array(name1, name2, name3, name4)

        do ival1 = 1, len(array1)
            call get_value(array1, ival1, array2, stat=stat)
            if (.not. allocated(ret)) then
                allocate (ret(len(array2), len(array1)))
            end if

            if (stat /= toml_stat%success) then
                call error_stop('String array 2d '//name1//'.'//name2//' is not found')
            end if

            do ival2 = 1, len(array2)
                call get_value(array2, ival1, ret(ival1, ival2)%string, stat=stat)

                if (stat /= toml_stat%success) then
                    call error_stop('String array 2d '//name1//'.'//name2//' is not found')
                end if
            end do
        end do
    end function

    function toml_to_string(self) result(ret)
        class(t_TomlWrapper), intent(inout) :: self
        character(len=:), allocatable :: ret

        ret = toml_serialize(self%table)
    end function

end module
