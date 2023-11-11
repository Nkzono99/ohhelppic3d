program main
    use M_CLI2, only: set_args, sget, set_mode, unnamed
    use ohhelppic3d

    implicit none

    character(:), allocatable :: parameter_filepath

    call parse_args()

    if (allocated(unnamed)) then
        if (size(unnamed) > 0) then
            parameter_filepath = unnamed(1)
        end if
    end if

    call pic(parameter_filepath)

contains

    subroutine parse_args
        call set_mode('strict')

        call set_args('')

        parameter_filepath = 'parameters.toml'
    end subroutine

end program main
