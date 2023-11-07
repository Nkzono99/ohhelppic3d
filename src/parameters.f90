module m_parameters
    use m_science_constants
    use m_toml_wrapper, only: t_TomlWrapper
    use m_string_holder, only: t_StringHolder
    implicit none

    private
    public t_Parameters, new_Parameters

    type :: t_Parameters
        character(len=:), allocatable :: toml_filepath
        type(t_TomlWrapper) :: toml

        character(len=:), allocatable :: simulation_name

        ! [simulation_type]
        character(len=:), allocatable :: simulation_type

        ! [continuaus]
        integer :: load_continuaus_data
        character(len=:), allocatable :: load_directory
        integer :: save_continuaus_data
        character(len=:), allocatable :: save_directory

        ! [system]
        integer :: nx
        integer :: ny
        integer :: nz
        integer :: nstep
        double precision :: dt
        integer :: nmacro_particles_per_grid(3)
        double precision :: particle_buffer_size(3)
        integer :: imbalance_tolerance_percentage = 10

        ! [system.method]
        character(len=:), allocatable :: particle_mover_type

        ! [system.mpi]
        integer :: nnodes(3)

        ! [system.outer_boundaries]
        type(t_StringHolder), allocatable :: boundary_communication(:)
        type(t_StringHolder), allocatable :: boundary_type_for_electromagnetic_field(:)
        type(t_StringHolder), allocatable :: boundary_type_for_particle(:, :)

        ! [output]
        ! [output.stdout]
        integer :: stdout_interval_step = 10

        ! [output.snapshot]
        integer :: output_start_step = 0
        integer :: field_output_interval
        integer :: current_output_interval

        ! [plasma]
        integer :: nspecies
        double precision, allocatable :: charge_to_mass_ratio(:)
        double precision, allocatable :: plasma_frequency(:)
        double precision, allocatable :: thermal_velocity_para(:)
        double precision, allocatable :: thermal_velocity_perp(:)
        double precision, allocatable :: flow_velocity(:)
        double precision, allocatable :: flow_angle_deg_z(:)
        double precision, allocatable :: flow_angle_deg_xy(:)
        double precision ::  cyclotron_frequency = 0.0

        type(t_StringHolder), allocatable :: plasma_initialization(:)

        double precision, allocatable :: charge_per_macro_particle(:)

    contains

    end type

contains

    function new_Parameters(toml_filepath) result(obj)
        character(len=*), intent(in) :: toml_filepath
        type(t_Parameters) :: obj

        obj%toml_filepath = toml_filepath
        call obj%toml%load(toml_filepath)

        ! TODO: Retrieve basic parameters.
        ! Parameters for complex settings are taken out in each module.
        obj%simulation_name = obj%toml%require_string1('simulation_name', default='main')

        ! [simulation_type]
        obj%simulation_type = obj%toml%require_string2('simulation_type', 'simulation_type')

        ! [continuaus]
        obj%load_continuaus_data = obj%toml%require_int('continuaus', 'load_continuaus_data')
        obj%load_directory = obj%toml%require_string2('continuaus', 'load_directory', default='SNAPSHOT0')
        obj%save_continuaus_data = obj%toml%require_int('continuaus', 'save_continuaus_data')
        obj%save_directory = obj%toml%require_string2('continuaus', 'save_directory', default='SNAPSHOT1')

        ! [system]
        obj%nx = obj%toml%require_int('system', 'nx')
        obj%ny = obj%toml%require_int('system', 'ny')
        obj%nz = obj%toml%require_int('system', 'nz')
        obj%nstep = obj%toml%require_int('system', 'nstep')
        obj%dt = obj%toml%require_double('system', 'dt')
        obj%nmacro_particles_per_grid = obj%toml%require_int_array('system', 'nmacro_particles_per_grid')
        obj%particle_buffer_size = obj%toml%require_double_array('system', 'particle_buffer_size')
        obj%imbalance_tolerance_percentage = obj%toml%require_int('system', 'imbalance_tolerance_percentage', default=10)

        ! [system.method]
        obj%particle_mover_type = obj%toml%require_string3('system', 'method', 'particle_mover', default='boris')

        ! [system.mpi]
        obj%nnodes = obj%toml%require_int_array('system', 'mpi', 'nnodes')

        ! [system.outer_boundaries]
        obj%boundary_communication = obj%toml%require_string_array('system', 'outer_boundary', 'boundary_communication')
        obj%boundary_type_for_electromagnetic_field = obj%toml%require_string_array('system', 'outer_boundary', 'boundary_type_for_electromagnetic_field')
        obj%boundary_type_for_particle = obj%toml%require_string_array2d('system', 'outer_boundary', 'boundary_type_for_particle')
        
        ! [output]
        ! [output.stdout]
        obj%stdout_interval_step = obj%toml%require_int('output', 'stdout', 'stdout_interval_step', default=10)

        ! [output.snapshot]
        obj%output_start_step = obj%toml%require_int('output', 'snapshot', 'output_start_step', default=0)
        obj%field_output_interval = obj%toml%require_int('output', 'snapshot', 'field_output_interval')
        obj%current_output_interval = obj%toml%require_int('output', 'snapshot', 'current_output_interval')

        ! [plasma]
        obj%nspecies = obj%toml%require_int('plasma', 'nspecies')
        obj%charge_to_mass_ratio = obj%toml%require_double_array('plasma', 'charge_to_mass_ratio')
        obj%plasma_frequency = obj%toml%require_double_array('plasma', 'plasma_frequency')
        obj%thermal_velocity_para = obj%toml%require_double_array('plasma', 'thermal_velocity_para')
        obj%thermal_velocity_perp = obj%toml%require_double_array('plasma', 'thermal_velocity_perp')
        obj%flow_velocity = obj%toml%require_double_array('plasma', 'flow_velocity')
        obj%flow_angle_deg_z = obj%toml%require_double_array('plasma', 'flow_angle_deg_z')
        obj%flow_angle_deg_xy = obj%toml%require_double_array('plasma', 'flow_angle_deg_xy')
        obj%cyclotron_frequency = obj%toml%require_double('plasma', 'cyclotron_frequency', default=0d0)

        ! [plasma.initialization]
        obj%plasma_initialization = obj%toml%require_string_array('plasma', 'initialization', 'plasma_initialization')

        block
            double precision, allocatable :: wp(:)
            double precision, allocatable :: qm(:)

            wp = obj%plasma_frequency(:)
            qm = obj%charge_to_mass_ratio(:)

            obj%charge_per_macro_particle = &
                wp*wp/qm/obj%nmacro_particles_per_grid
        end block
    end function

end module
