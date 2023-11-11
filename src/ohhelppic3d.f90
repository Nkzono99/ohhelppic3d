module ohhelppic3d
    use mpi
    use m_ohhelp, only: t_OhHelp, new_OhHelp
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    use m_ohfield_factory, only: t_OhFieldFactory, new_OhFieldFactory
    use m_ohfield, only: t_OhField, tp_OhField, &
                         BOUNDARY_CONDITION_PERIODIC, &
                         BOUNDARY_CONDITION_NO_PERIODIC
    use m_ohparticles, only: t_OhParticles, new_OhParticles
    use oh_type, only: oh_particle
    use m_hdf5
    use m_particle_mover, only: t_ParticleMover
    use m_particle_mover_factory, only: t_ParticleMoverFactory, new_ParticleMoverFactory
    use m_parameters, only: t_Parameters, new_Parameters
    use m_particle_injector_manager, only: t_ParticleInjectorManager, new_ParticleInjectorManager
    use m_no_position_distribution, only: new_NoPositionDistribution3d
    use m_velocity_distribution, only: new_NoVelocityDistribution3d
    use m_interpolator, only: t_Interpolator
    use m_mpi_fft_solver, only: t_MPIFFTSolver3d
    use m_mpi_fft_solver_factory
    use m_block, only: t_Block, new_Block
    use m_random_generator
    use m_pcg_generator
    use m_field_solver
    use m_poisson_field_solver
    use m_field_boundary_type
    use m_linear_interpolator
    use m_particle_boundaries
    use m_scatter
    use m_linear_scatter
    use m_str
    use m_hdf5_for_ohfield
    use m_domain
    implicit none

    private
    public pic

    !! Main Variables
    integer :: nprocs, myid
    type(t_OhHelp) :: ohhelp
    type(t_OhParticles) :: ohparticles
    integer :: pbase(3)
    type(t_OhField), target :: eb, aj, rho, phi
    class(t_ParticleInjectorManager), allocatable :: particle_injector_manager
    class(t_ParticleMover), allocatable :: particle_mover
    class(t_MPIFFTSolver3d), target, allocatable :: mpifft_solver3d
    class(t_Interpolator), allocatable :: interpolator
    type(t_Parameters) :: parameters
    class(t_RandomGenerator), allocatable, target :: random_generator
    character(len=15) :: toml_filepath = 'parameters.toml'
    class(t_FieldSolver), allocatable, target :: field_solver
    class(t_ParticleBoundaries), allocatable :: particle_boundaries
    class(t_Scatter), allocatable :: scatter
    class(t_Hdf5ForOhfield), allocatable :: hdf5_phi
    class(t_Hdf5ForOhfield), allocatable :: hdf5_rho

contains

    subroutine pic(parameter_filepath)
        character(*), intent(in) :: parameter_filepath
        integer :: istep

        block
            integer :: ierr

            call mpi_init(ierr)
            if (ierr /= 0) error stop "mpi_init failed"

            call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
            call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
            if (ierr /= 0) error stop "mpi_comm_rank failed"

            print *, myid, '/', nprocs
        end block

        parameters = new_Parameters(parameter_filepath)

        call initialize
        call initialize_model

        block
            integer :: ierr
            call hdf5_initialize(ierr)
        end block

        hdf5_phi = new_Hdf5ForOhfield('phisp00_0000.h5', 'phisp', &
                                      int([parameters%nx + 1, parameters%ny + 1, parameters%nz + 1], kind=8), &
                                      int([0, 0, 0], kind=8), &
                                      MPI_COMM_WORLD)
        hdf5_rho = new_Hdf5ForOhfield('rho00_0000.h5', 'rho', &
                                      int([parameters%nx + 1, parameters%ny + 1, parameters%nz + 1], kind=8), &
                                      int([0, 0, 0], kind=8), &
                                      MPI_COMM_WORLD)

        call hdf5_phi%write('0', phi, 1)
        call hdf5_rho%write('0', rho, 1)

        do istep = 1, parameters%nstep
            if (myid == 0 .and. mod(istep, parameters%stdout_interval_step) == 0) then
                print *, '------ '//str(istep)//' -------'
            end if

            call mainstep(parameters%dt)

            if (mod(istep, parameters%field_output_interval) == 0) then
                call hdf5_phi%write(str(istep/parameters%field_output_interval), phi, 1)
                call hdf5_rho%write(str(istep/parameters%field_output_interval), rho, 1)
            end if

            block
                integer :: ierr
                call MPI_Barrier(MPI_COMM_WORLD, ierr)
            end block
        end do

        call hdf5_phi%close()
        call hdf5_rho%close()

        block
            integer :: ierr
            call MPI_Barrier(MPI_COMM_WORLD, ierr)

            call hdf5_finalize(ierr)
            call mpi_finalize(ierr)
        end block

        call finalize
    end subroutine

    subroutine initialize
        random_generator = new_PcgGenerator([int(42, kind=8), int(52, kind=8)])
        call random_generator%advance(myid*int(100000000, kind=8))

        ! Init ohfields
        block
            type(t_OhFieldFactory) :: ohfield_factory
            ohfield_factory = new_OhFieldFactory()

            eb = ohfield_factory%create_field('electromagnetic')
            aj = ohfield_factory%create_field('current')
            rho = ohfield_factory%create_field('density')
            phi = ohfield_factory%create_field('potential')
        end block

        ! Init ohparticles
        block
            integer(kind=8) :: pbuf_size(parameters%nspecies)
            integer :: nx, ny, nz

            nx = parameters%nx; ny = parameters%ny; nz = parameters%nz
            pbuf_size(:) = parameters%nmacro_particles_per_grid(:) &
                           *parameters%particle_buffer_size(:) &
                           *nx*ny*nz
            ohparticles = new_OhParticles(parameters%nspecies, sum(pbuf_size), product(parameters%nnodes))
        end block

        ! new ohhelp
        block
            integer :: boundary_conditions(2, 3)
            integer :: i
            integer :: ierr
            character(len=:), allocatable :: name

            do i = 1, 3
                name = parameters%boundary_communication(i)%string
                if (name == 'periodic') then
                    boundary_conditions(:, i) = BOUNDARY_CONDITION_PERIODIC
                else if (name == 'no_periodic') then
                    boundary_conditions(:, i) = BOUNDARY_CONDITION_NO_PERIODIC
                else
                    write (stderr, '(a)') 'Error: '//'system.outer_boundary.boundary_communication is invalid: '//name
                    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
                    stop 1
                end if
            end do

            ohhelp = new_OhHelp(parameters%nspecies, &
                                parameters%nnodes, &
                                parameters%nx, parameters%ny, parameters%nz, &
                                boundary_conditions, &
                                parameters%imbalance_tolerance_percentage)
        end block

        ! Init ohhelp
        block
            type(tp_OhField) :: ohfields(4)

            ohfields(1)%ref => eb
            ohfields(2)%ref => aj
            ohfields(3)%ref => rho
            ohfields(4)%ref => phi

            call ohhelp%initialize(ohparticles, ohfields)
        end block

        ! Init particle injector
        particle_injector_manager = new_ParticleInjectorManager(parameters, random_generator)
        scatter = new_LinearScatter()

        block
            type(t_Block) :: local_block
            type(t_Block) :: global_block
            double precision :: boundary_values(2, 3)

            local_block = new_Block(ohhelp%subdomain_range(1, :, ohhelp%subdomain_id(1) + 1), &
                                    ohhelp%subdomain_range(2, :, ohhelp%subdomain_id(1) + 1))
            global_block = new_Block([0, 0, 0], [parameters%nx, parameters%ny, parameters%nz])
            boundary_values = reshape([[0d0, 0d0], [0d0, 0d0], [0d0, 0d0]], [2, 3])

            field_solver = new_PoissonFieldSolver(local_block, global_block, &
                                                  'fftw3', &
                                                  parameters%boundary_type_for_electromagnetic_field, &
                                                  boundary_values, &
                                                  myid, &
                                                  nprocs, &
                                                  MPI_COMM_WORLD, tag=10)
        end block

        interpolator = new_LinearInterpolator()

        block
            type(t_ParticleMoverFactory) :: particle_mover_factory

            particle_mover_factory = new_ParticleMoverFactory()
            particle_mover = particle_mover_factory%create_particle_mover(parameters%particle_mover_type)
        end block

        particle_boundaries = new_ParticleBoundaries(parameters%nx, parameters%ny, parameters%nz)

    end subroutine

    subroutine initialize_model()
        call particle_injector_manager%initialize_particles(ohhelp)

        ! Correct load balancing
        block
            type(tp_OhField) :: ohfields(4)

            ohfields(1)%ref => eb
            ohfields(2)%ref => aj
            ohfields(3)%ref => rho
            ohfields(4)%ref => phi

            call ohhelp%correct_load_balancing(ohparticles, eb, ohfields_to_be_notified=ohfields)
        end block

        call scatter_charge

        call field_solver%solve(0d0, rho, aj, eb, phi, ohhelp)
    end subroutine

    subroutine mainstep(dt)
        double precision, intent(in) :: dt

        call update_particle(dt)

        call scatter_charge

        call field_solver%solve(dt, rho, aj, eb, phi, ohhelp)

        call particle_injector_manager%inject_particles(dt, ohhelp)

        block
            type(tp_OhField) :: ohfields(4)

            ohfields(1)%ref => eb
            ohfields(2)%ref => aj
            ohfields(3)%ref => rho
            ohfields(4)%ref => phi

            call ohhelp%correct_load_balancing(ohparticles, eb, ohfields_to_be_notified=ohfields)
        end block

    end subroutine

    subroutine scatter_charge
        rho%values(:, :, :, :, :) = 0

        block
            integer :: eps
            integer :: ps
            integer :: ispec, ipcl
            integer :: ipcl_start, ipcl_end

            if (ohhelp%is_primary_mode()) then
                eps = 1
            else
                eps = 2
            end if

            do ps = 1, eps
                do ispec = 1, ohparticles%nspecies

                    ipcl_start = ohparticles%start_index(ispec, ps)
                    ipcl_end = ohparticles%end_index(ispec, ps)
                    do ipcl = ipcl_start, ipcl_end
                        call scatter%scatter(ohparticles%pbuf(ipcl), rho, [parameters%charge_per_macro_particle(ispec)], ps)
                    end do
                end do
            end do
        end block

        if (ohhelp%is_secondary_mode()) then
            call ohhelp%reduce_field(rho)
        end if

        call ohhelp%exchange_borders(rho)

        block
            integer :: eps
            integer :: ps

            if (ohhelp%is_primary_mode()) then
                eps = 1
            else
                eps = 2
            end if

            do ps = 1, eps
                block
                    integer :: xl, yl, zl
                    integer :: xu, yu, zu

                    xl = 0
                    yl = 0
                    zl = 0
                    xu = rho%subdomain_range(2, 1, ps) - rho%subdomain_range(1, 1, ps)
                    yu = rho%subdomain_range(2, 2, ps) - rho%subdomain_range(1, 2, ps)
                    zu = rho%subdomain_range(2, 3, ps) - rho%subdomain_range(1, 3, ps)

                    rho%values(1, xl, :, :, ps) = rho%values(1, xl, :, :, ps) + rho%values(1, xl - 1, :, :, ps)
                    rho%values(1, :, yl, :, ps) = rho%values(1, :, yl, :, ps) + rho%values(1, :, yl - 1, :, ps)
                    rho%values(1, :, :, zl, ps) = rho%values(1, :, :, zl, ps) + rho%values(1, :, :, zl - 1, ps)
                    rho%values(1, xu, :, :, ps) = rho%values(1, xu, :, :, ps) + rho%values(1, xu + 1, :, :, ps)
                    rho%values(1, :, yu, :, ps) = rho%values(1, :, yu, :, ps) + rho%values(1, :, yu + 1, :, ps)
                    rho%values(1, :, :, zu, ps) = rho%values(1, :, :, zu, ps) + rho%values(1, :, :, zu + 1, ps)
                end block
            end do
        end block
    end subroutine

    subroutine update_particle(dt)
        double precision, intent(in) :: dt

        integer :: ps
        integer :: ispec

        integer :: eps

        if (ohhelp%is_primary_mode()) then
            eps = 1
        else
            eps = 2
        end if

        do ps = 1, eps
        do ispec = 1, ohparticles%nspecies
            block
                double precision :: eb_interped(6) = 0
                integer :: ipcl_start, ipcl_end
                integer :: ipcl
                double precision :: qm

                qm = parameters%charge_to_mass_ratio(ispec)
                ipcl_start = ohparticles%start_index(ispec, ps)
                ipcl_end = ohparticles%end_index(ispec, ps)

                do ipcl = ipcl_start, ipcl_end
                    eb_interped(:) = interpolator%interp(ohparticles%pbuf(ipcl), eb, ps)

                    call particle_mover%move(ohparticles%pbuf(ipcl), qm, eb_interped, dt)

                    call particle_boundaries%apply(ohparticles%pbuf(ipcl), dt)

                    call ohhelp%correct_particle(ohparticles%pbuf(ipcl), ps)
                end do
            end block
        end do
        end do
    end subroutine

    subroutine finalize
        integer :: status
    end subroutine

end module ohhelppic3d
