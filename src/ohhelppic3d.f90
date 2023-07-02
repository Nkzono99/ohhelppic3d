module ohhelppic3d
    use mpi
    use m_ohhelp, only: t_OhHelp, new_OhHelp
    use m_ohfield_factory, only: t_OhFieldFactory, new_OhFieldFactory
    use m_ohfield, only: t_OhField, tp_OhField, &
                         BC_PERIODIC => BOUNDARY_CONDITION_PERIODIC, &
                         BC_NO_PERIODIC => BOUNDARY_CONDITION_NO_PERIODIC
    use m_ohparticles, only: t_OhParticles, new_OhParticles
    use oh_type, only: oh_particle
    use m_hdf5
    use m_particle_mover, only: t_ParticleMover
    use m_particle_mover_factory, only: t_ParticleMoverFactory, new_ParticleMoverFactory
    implicit none
    private
    public pic

    !! Main Variables
    integer :: nprocs, myid
    type(t_OhHelp) :: ohhelp
    type(t_OhParticles) :: ohparticles
    integer :: pbase(3)
    type(t_OhField), target :: eb, aj, rho, phi
    class(t_ParticleMover), allocatable :: particle_mover

    !! Parameters
    integer :: nspec = 2
    integer :: tolerance = 10
    integer(kind=8) :: max_npcls = 1000
    integer :: nnodes(3) = [2, 1, 2]

    integer :: boundary_conditions(2, 3) = [[BC_PERIODIC, BC_PERIODIC], &
                                            [BC_PERIODIC, BC_PERIODIC], &
                                            [BC_PERIODIC, BC_PERIODIC]]

    integer :: nx = 8
    integer :: ny = 4
    integer :: nz = 16

    integer :: nstep = 1000
    double precision :: dt = 0.1d0

contains

    subroutine pic

        integer :: istep

        ! TODO:パラメータの読み込みの実装

        call initialize

        do istep = 1, nstep
            call mainstep(dt)
        end do

        call finalize
    end subroutine

    subroutine initialize
        integer :: status
        integer :: ierr
        type(t_ParticleMoverFactory) :: particle_mover_factory

        call mpi_init(ierr)
        if (ierr /= 0) error stop "mpi_init failed"

        call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
        call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
        if (ierr /= 0) error stop "mpi_comm_rank failed"
        print *, myid, '/', nprocs

        call hdf5_initialize(status)

        call ohinit

        particle_mover_factory = new_ParticleMoverFactory()
        call particle_mover_factory%create_particle_mover('boris', particle_mover)

        ! TODO: 単位系のリスケール

        ! TODO: 粒子の初期化の実装

        ! TODO: 場の初期化
    end subroutine

    subroutine mainstep(dt)
        integer :: status
        double precision, intent(in) :: dt

        ! TODO: 粒子の生成

        ! TODO: 場の更新 & 場のリロケート(self-forceの回避)

        ! TODO: 粒子の更新

        ! TODO: 粒子密度 or 電流の場への配分

        ! TODO: ロードバランスの適用

        ! TODO: スナップショットの出力

    end subroutine

    subroutine finalize
        integer :: status

        call hdf5_finalize(status)
    end subroutine

    ! subroutine ohinit
    !     type(t_OhFieldFactory) :: ohfield_factory

    !     integer :: ii(5)

    !     ohfield_factory = new_OhFieldFactory()

    !     eb = ohfield_factory%create_field('electromagnetic', 2)
    !     aj = ohfield_factory%create_field('current', 2)
    !     rho = ohfield_factory%create_field('density', 2)
    !     phi = ohfield_factory%create_field('potential', 2)

    !     ohparticles = new_OhParticles(nspec, max_npcls, product(nnodes), tolerance)

    !     ohhelp = new_OhHelp(nnodes, nx, ny, nz, boundary_conditions)

    !     call ohhelp%initialize(ohparticles, [eb, aj, rho, phi])

    !     call ohhelp%transbound(ohparticles)
    ! end subroutine

    subroutine ohinit
        type(t_OhFieldFactory) :: ohfield_factory

        integer :: ii(5)
        integer :: i
        type(oh_particle) :: particle
        integer :: subdomain(2, 3)

        ohfield_factory = new_OhFieldFactory()

        eb = ohfield_factory%create_field('electromagnetic', 2)
        aj = ohfield_factory%create_field('current', 2)
        rho = ohfield_factory%create_field('density', 2)
        phi = ohfield_factory%create_field('potential', 2)

        ohparticles = new_OhParticles(nspec, max_npcls, product(nnodes))

        ohhelp = new_OhHelp(nnodes, nx, ny, nz, boundary_conditions, tolerance)

        block ! Initialization of OhHelp
            type(tp_OhField) :: ohfields(4)
            ohfields(1)%ref => eb
            ohfields(2)%ref => aj
            ohfields(3)%ref => rho
            ohfields(4)%ref => phi

            call ohhelp%initialize(ohparticles, ohfields)
            call ohhelp%transbound(ohparticles)
        end block

        do i = 1, 30*(ohhelp%subdomain_id(1)+1)
            subdomain = ohhelp%subdomain_range(:, :, ohhelp%subdomain_id(1)+1)
            block
                double precision :: position(3)

                position(:) = subdomain(1, :) + (subdomain(2, :) - subdomain(1, :))/1000d0*i
                particle%x = position(1)
                particle%y = position(2)
                particle%z = position(3)
                particle%nid = ohhelp%subdomain_id(1)
                particle%pid = 0
                particle%spec = 1
            end block

            call ohhelp%inject_particle(particle)
        end do

        call ohhelp%transbound(ohparticles)
    end subroutine

end module ohhelppic3d
