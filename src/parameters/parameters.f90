module m_parameters
    use m_toml_wrapper, only: t_TomlWrapper
    implicit none

    private
    public t_Parameters, new_Parameters

    type :: t_Parameters
        character(len=:), allocatable :: toml_filepath
        type(t_TomlWrapper) :: toml

        character(len=:), allocatable :: simulation_name

        integer :: nx
        integer :: ny
        integer :: nz

        integer :: nstep
        double precision :: dt

        integer :: nspecies

        integer :: boundary_conditions(2, 3)

        integer, allocatable :: max_npcls(:) ! (nspecies)

        integer :: tolerance
        integer :: nnodes(3)

    contains

    end type

contains

    function new_Parameters(toml_filepath) result(obj)
        character(len=*), intent(in) :: toml_filepath
        type(t_Parameters) :: obj

        obj%toml_filepath = toml_filepath
        call obj%toml%load(toml_filepath)

        ! Retrieve basic parameters.
        ! Parameters for complex settings are taken out in each module.
        obj%simulation_name = obj%toml%require_string1('simulation_name')
        obj%nx = obj%toml%require_int('system', 'nx')
        obj%ny = obj%toml%require_int('system', 'ny')
        obj%nz = obj%toml%require_int('system', 'nz')
        obj%nstep = obj%toml%require_int('system', 'nstep')
        obj%dt = obj%toml%require_double('system', 'dt')
        obj%nspecies = obj%toml%require_int('plasma', 'nspecies')
    end function

end module
