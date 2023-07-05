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
    use m_position_distribution, only: new_NoPositionDistribution3d
    use m_velocity_distribution, only: new_NoVelocityDistribution3d
    use m_interpolator, only: t_Interpolator
    use m_mpi_fft_solver, only: t_MPIFFTSolver3d
    use m_mpi_fftw_solver, only: new_MPIFFTWSolver3d
    use m_block, only: t_Block, new_Block
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
    class(t_MPIFFTSolver3d), allocatable :: mpifft_solver3d
    class(t_Interpolator), allocatable :: interpolator
    type(t_Parameters) :: parameters
    character(len=15) :: toml_filepath = 'parameters.toml'

contains

    subroutine pic
        integer :: istep
        integer :: ierr

        call mpi_init(ierr)
        if (ierr /= 0) error stop "mpi_init failed"

        call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
        call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
        if (ierr /= 0) error stop "mpi_comm_rank failed"
        print *, myid, '/', nprocs

        call initialize

        do istep = 1, parameters%nstep
            call mainstep(parameters%dt)
        end do

        call finalize
    end subroutine

    subroutine initialize
        integer :: status

        call hdf5_initialize(status)

        parameters = new_Parameters(toml_filepath)

        block
            type(t_OhFieldFactory) :: ohfield_factory
            ohfield_factory = new_OhFieldFactory()

            eb = ohfield_factory%create_field('electromagnetic')
            aj = ohfield_factory%create_field('current')
            rho = ohfield_factory%create_field('density')
            phi = ohfield_factory%create_field('potential')
        end block

        block
            integer(kind=8) :: pbuf_size(parameters%nspecies)
            integer :: nx, ny, nz

            nx = parameters%nx; ny = parameters%ny; nz = parameters%nz
            pbuf_size(:) = parameters%nmacro_particles_per_grid(:) &
                           *parameters%particle_buffer_size(:) &
                           *nx*ny*nz
            ohparticles = new_OhParticles(parameters%nspecies, sum(pbuf_size), product(parameters%nnodes))
        end block

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

            ohhelp = new_OhHelp(parameters%nnodes, &
                                parameters%nx, parameters%ny, parameters%nz, &
                                boundary_conditions, &
                                parameters%imbalance_tolerance_percentage)
        end block

        block
            type(tp_OhField) :: ohfields(4)

            ohfields(1)%ref => eb
            ohfields(2)%ref => aj
            ohfields(3)%ref => rho
            ohfields(4)%ref => phi

            call ohhelp%initialize(ohparticles, ohfields)
        end block

        ! TODO: 粒子の初期化の実装
        particle_injector_manager = new_ParticleInjectorManager()
        call particle_injector_manager%initialize_particles(ohhelp, parameters)

        block
            type(tp_OhField) :: ohfields(4)

            ohfields(1)%ref => eb
            ohfields(2)%ref => aj
            ohfields(3)%ref => rho
            ohfields(4)%ref => phi

            call ohhelp%correct_load_balancing(ohparticles, eb, ohfields_to_be_notified=ohfields)
        end block

        ! TODO: 場の初期化
        ! TODO: 粒子Scatterの実装

        ! TODO: フィールドソルバーの実装
        block
            integer :: fft_boundary_types(3)
            type(t_Block) :: local_block

            fft_boundary_types = [0, 0, 0]

            local_block = new_Block(ohhelp%subdomain_range(1, :, ohhelp%subdomain_id(1)) + 1, &
                                    ohhelp%subdomain_range(2, :, ohhelp%subdomain_id(1)) + 1)

            mpifft_solver3d = new_MPIFFTWSolver3d(fft_boundary_types, &
                                                  local_block, &
                                                  parameters%nx, parameters%ny, parameters%nz, &
                                                  myid + 1, nprocs, &
                                                  MPI_COMM_WORLD, tag=10)
        end block

        block
            integer :: lnx, lny, lnz

            lnx = rho%subdomain_range(2, 1, 1) - rho%subdomain_range(1, 1, 1) + 1
            lnx = rho%subdomain_range(2, 2, 1) - rho%subdomain_range(1, 2, 1) + 1
            lnx = rho%subdomain_range(2, 3, 1) - rho%subdomain_range(1, 3, 1) + 1

            rho%values(1, 2, 2, 2, 1) = 2.125d0

            call mpifft_solver3d%forward(rho%values(1, 0:lnx, 0:lnx, 0:lnx, 1), &
                                         rho%values(1, 0:lnx, 0:lnx, 0:lnx, 1))
            call mpifft_solver3d%backward(rho%values(1, 0:lnx, 0:lnx, 0:lnx, 1), &
                                         rho%values(1, 0:lnx, 0:lnx, 0:lnx, 1))

            print *, rho%values(1, 2, 2, 2, 1)
        end block

        block
            type(t_ParticleMoverFactory) :: particle_mover_factory

            particle_mover_factory = new_ParticleMoverFactory()
            particle_mover = particle_mover_factory%create_particle_mover(parameters%particle_mover_type)
        end block
    end subroutine

    subroutine mainstep(dt)
        double precision, intent(in) :: dt

        call particle_injector_manager%inject_particles(dt, ohhelp)

        ! TODO: 場の更新 & 場のリロケート(self-forceの回避)

        ! TODO: 粒子の更新
        block
            integer :: ps
            integer :: ispec, ipcl
            integer :: ipcl_start, ipcl_end
            double precision :: qm
            double precision :: eb_interped(6)

            ! Primaryモードの場合は、start_index(ispec, 2) > end_index(ispec, 2)のため最内ループは実行されない
            do ps = 1, 2
            do ispec = 1, ohparticles%nspecies
                ipcl_start = ohparticles%start_index(ispec, ps)
                ipcl_end = ohparticles%end_index(ispec, ps)
                do ipcl = ipcl_start, ipcl_end
                    eb_interped(:) = interpolator%interp(ohparticles%pbuf(ipcl), eb)

                    call particle_mover%move(ohparticles%pbuf(ipcl), qm, eb_interped, dt)
                end do
            end do
            end do
        end block

        ! TODO: 粒子密度 or 電流の場への配分

        block
            type(tp_OhField) :: ohfields(4)

            ohfields(1)%ref => eb
            ohfields(2)%ref => aj
            ohfields(3)%ref => rho
            ohfields(4)%ref => phi

            call ohhelp%correct_load_balancing(ohparticles, eb, ohfields_to_be_notified=ohfields)
        end block

        ! TODO: スナップショットの出力

    end subroutine

    subroutine finalize
        integer :: status

        call hdf5_finalize(status)
    end subroutine

end module ohhelppic3d
