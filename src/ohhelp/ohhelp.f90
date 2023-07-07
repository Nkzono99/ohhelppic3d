#define OH_LIB_LEVEL 3

#ifndef OH_PBUF_SIZE
#define OH_PBUF_SIZE 16384
#endif

!> OhHelp Wrapper module.
module m_ohhelp
    use ohhelp2, only: oh2_max_local_particles, oh2_inject_particle
    use ohhelp3, only: oh3_init, oh3_transbound, &
                       oh3_map_particle_to_neighbor, &
                       oh3_map_particle_to_subdomain, &
                       oh3_bcast_field
    use oh_type, only: oh_particle, oh_mycomm
    use m_ohparticles, only: t_OhParticles
    use m_ohfield, only: t_FieldExtensionInfo, new_FieldExtensionInfo, &
                         tp_OhField, &
                         t_BoundaryCommunicationInfo, new_BoundaryCommunicationInfo, &
                         t_BoundaryCommunicationInfos, new_BoundaryCommunicationInfos, &
                         t_OhField, new_OhField, &
                         NBOUNDARY_CONDITION_TYPES, &
                         BOUNDARY_CONDITION_PERIODIC, BOUNDARY_CONDITION_NO_PERIODIC
    implicit none

    integer, parameter :: PARTICLE_BUFSIZE = OH_PBUF_SIZE

    !> Helpand-helper configuration in (re)build and in secondary mode (= -1).
    integer, parameter :: OHH_REQUIRES_BCAST = -1
    !> In primary mode (= 0).
    integer, parameter :: OHH_PRIMARY_MODE = 0
    !> In secondary mode (= 1).
    integer, parameter :: OHH_SECONDARY_MODE = 1

    private
    public oh_particle
    public t_OhHelp, new_OhHelp
    ! integer(kind=4), external :: oh2_max_local_particles
    ! integer(kind=4), external :: oh3_transbound
    ! integer(kind=4), external :: oh3_map_particle_to_neighbor
    ! integer(kind=4), external :: oh3_map_particle_to_subdomain

    type t_OhHelp
        integer :: subdomain_id(2)
        !> OHH_PRIMARY_MODE: primary mode, OHH_SECONDARY_MODE: secondary mode, OHH_REQUIRES_BCAST: requires bcast
        integer :: current_mode = 0

        integer :: loadbalance_tolerance_percentage = 10

        type(oh_mycomm) :: communicator

        integer :: neighber_subdomain_ids(3, 3, 3)
        integer :: process_coordinates(3)
        integer, allocatable :: subdomain_range(:, :, :) ! (2, 3, nprocs)
        integer :: whole_domain_range(2, 3)

        integer :: nboundary_condition_types = NBOUNDARY_CONDITION_TYPES
        integer :: boundary_conditions(2, 3) = 1
        integer, allocatable :: subdomain_boundary_conditions(:, :, :) ! (nprocs, 3, 2)

        type(t_FieldExtensionInfo), allocatable :: field_extension_infos(:) ! (nextension + 1)
        type(t_BoundaryCommunicationInfos), allocatable :: boundary_communication_infos(:)

        integer, allocatable :: field_sizes(:, :, :) ! (2, 3, F)

    contains

        procedure :: initialize => ohhelp_initialize
        procedure :: allocate_ohfield => ohhelp_allocate_ohfield

        procedure :: correct_load_balancing => ohhelp_correct_load_balancing

        procedure :: transbound => ohhelp_transbound
        procedure :: inject_particle => ohhelp_inject_particle

        procedure :: requires_broadcast_field => ohhelp_requires_broadcast_field
        procedure :: broadcast_field => ohhelp_broadcast_field

        procedure :: is_primary_mode => ohhelp_is_primary_mode
        procedure :: is_secondary_mode => ohhelp_is_secondary_mode

        procedure, private :: set_field_extension_infos => ohhelp_set_field_extension_infos
        procedure, private :: set_boundary_communication_infos => ohhelp_set_boundary_communication_infos

        procedure :: notify_subdomain_id => ohhelp_notify_subdomain_id
    end type

contains

    function new_OhHelp(nnodes, &
                        nx, ny, nz, &
                        boundary_conditions, &
                        loadbalance_tolerance_percentage) result(obj)
        integer, intent(in) :: nnodes(3)
        integer, intent(in) :: nx
        integer, intent(in) :: ny
        integer, intent(in) :: nz
        integer, intent(in) :: loadbalance_tolerance_percentage

        !> Boundary conditions.
        !>
        !>   boundary_conditions :=
        !>     [[xlower, xupper], [ylower, yupper], [zlower, zupper]].
        !>
        !>   [x/y/z][lower/upper] :=
        !>     periodic boundary : BOUNDARY_CONDITION_PERIODIC
        !>         otherwise     : BOUNDARY_CONDITION_NO_PERIODIC)
        integer, intent(in) :: boundary_conditions(2, 3)
        type(t_OhHelp) :: obj

        integer :: nprocs

        nprocs = product(nnodes)

        allocate (obj%subdomain_range(2, 3, nprocs))
        allocate (obj%subdomain_boundary_conditions(2, 3, nprocs))

        obj%process_coordinates(:) = nnodes(:)
        obj%whole_domain_range(:, :) = reshape([0, nx, 0, ny, 0, nz], [2, 3])
        obj%boundary_conditions(:, :) = boundary_conditions(:, :)
        obj%loadbalance_tolerance_percentage = loadbalance_tolerance_percentage
    end function

    subroutine ohhelp_initialize(self, ohparticles, ohfields)
        class(t_OhHelp), intent(inout) :: self
        type(t_OhParticles), intent(inout) :: ohparticles
        type(tp_OhField), intent(inout) :: ohfields(:)

        integer :: status = 0
        integer :: repiter = 0
        integer :: verbose = 0

        integer, allocatable :: ftypes(:, :)
        integer, allocatable :: cfields(:)
        integer, allocatable :: ctypes(:, :, :, :)

        call ohparticles%allocate_pbuf(oh2_max_local_particles(ohparticles%max_nparticles, &
                                                               self%loadbalance_tolerance_percentage, &
                                                               PARTICLE_BUFSIZE))

        call self%set_field_extension_infos(ohfields)
        call self%set_boundary_communication_infos(ohfields)

        allocate (ftypes(7, size(self%field_extension_infos) + 1))
        allocate (cfields(size(self%boundary_communication_infos) + 1))
        allocate (ctypes(3, 2, NBOUNDARY_CONDITION_TYPES, size(self%boundary_communication_infos)))

        ftypes = ftypes_from(self%field_extension_infos)
        cfields = cfields_from(self%boundary_communication_infos)
        ctypes = ctypes_from(self%boundary_communication_infos)

        ! Let ohhelp set up the following parameters.
        self%neighber_subdomain_ids(1, 1, 1) = -1
        self%subdomain_range(1, 1, 1) = 0
        self%subdomain_range(2, 1, 1) = -1

        call oh3_init(self%subdomain_id(:), & ! sdid(2)
                      ohparticles%nspecies, & ! nspec
                      self%loadbalance_tolerance_percentage, & ! maxfrac
                      ohparticles%particle_count_histgram(:, :, :), & ! nphgram(N, S, 2)
                      ohparticles%total_local_particles(:, :), & ! totalp(S, 2)
                      ohparticles%pbuf(:), & ! pbuf(:)
                      ohparticles%pbase(:), & ! pbase(3)
                      ohparticles%max_local_particles, & ! maxlocalp
                      self%communicator, & ! mycomm
                      self%neighber_subdomain_ids(:, :, :), & ! nbor(3, 3, 3)
                      self%process_coordinates(:), & ! pcoord(3)
                      self%subdomain_range(:, :, :), & ! sdoms(2, 3, N)
                      self%whole_domain_range(:, :), & ! scoord(2, 3)
                      self%nboundary_condition_types, & ! nbound
                      self%boundary_conditions(:, :), & ! bcond(2, 3)
                      self%subdomain_boundary_conditions(:, :, :), & ! bounds(2, 3, N)
                      ftypes(:, :), & ! ftypes(7, F+:)
                      cfields(:), & ! cfields(C+:)
                      ctypes(:, :, :, :), & ! ctypes(3, 2, B, C)
                      self%field_sizes(:, :, :), & ! fsizes(2, 3, F)
                      status, &
                      repiter, &
                      verbose)

        ! Initialize particle-related variables, etc. in ohhelp (first time calling transbound function)
        call self%transbound(ohparticles)

        block
            integer :: iohfield
            do iohfield = 1, size(ohfields)
                call self%allocate_ohfield(ohfields(iohfield)%ref)
                call self%notify_subdomain_id(ohfields(iohfield)%ref)
            end do
        end block

    contains

        function ftypes_from(infos) result(ret)
            type(t_FieldExtensionInfo), intent(in) :: infos(:)
            integer :: ret(7, size(infos) + 1)

            block
                integer :: i, id

                do i = 1, size(infos)
                    id = infos(i)%id
                    ret(1, id) = infos(i)%nelements
                    ret(2:3, id) = infos(i)%nextensions(1:2)
                    ret(4:5, id) = infos(i)%nextensions_for_broadcast(1:2)
                    ret(6:7, id) = infos(i)%nextensions_for_reduction(1:2)
                end do
            end block

            ret(1, size(infos) + 1) = 0

        end function

        function cfields_from(infos) result(ret)
            type(t_BoundaryCommunicationInfos), intent(in) :: infos(:)
            integer :: ret(size(infos) + 1)

            ret(1:size(infos)) = infos(:)%id
            ret(size(infos) + 1) = 0
        end function

        function ctypes_from(infos) result(ret)
            type(t_BoundaryCommunicationInfos), intent(in) :: infos(:)
            integer :: ret(3, 2, NBOUNDARY_CONDITION_TYPES, size(infos))

            integer :: i, ib

            do i = 1, size(infos)
                do ib = 1, NBOUNDARY_CONDITION_TYPES
                    ret(1, 1, ib, i) = infos(i)%infos(ib)%downward_comm_send_offset
                    ret(2, 1, ib, i) = infos(i)%infos(ib)%downward_comm_receive_offset
                    ret(3, 1, ib, i) = infos(i)%infos(ib)%downward_comm_nsends
                    ret(1, 2, ib, i) = infos(i)%infos(ib)%upward_comm_send_offset
                    ret(2, 2, ib, i) = infos(i)%infos(ib)%upward_comm_receive_offset
                    ret(3, 2, ib, i) = infos(i)%infos(ib)%upward_comm_nsends
                end do
            end do
        end function

    end subroutine

    subroutine ohhelp_set_field_extension_infos(self, ohfields)
        class(t_OhHelp), intent(inout) :: self
        type(tp_OhField), intent(in) :: ohfields(:)

        allocate (self%field_extension_infos(size(ohfields)))
        allocate (self%field_sizes(2, 3, size(ohfields)))

        block
            integer :: i
            integer :: id
            type(t_FieldExtensionInfo) :: info
            do i = 1, size(ohfields)
                info = ohfields(i)%ref%extension_info
                id = info%id
                self%field_extension_infos(id) = info
            end do
        end block
    end subroutine

    subroutine ohhelp_set_boundary_communication_infos(self, ohfields)
        class(t_OhHelp), intent(inout) :: self
        type(tp_OhField), intent(in) :: ohfields(:)

        allocate (self%boundary_communication_infos(size(ohfields)))

        block
            integer :: i, j
            integer :: id
            type(t_BoundaryCommunicationInfos), allocatable :: infos(:)
            do i = 1, size(ohfields)
                infos = ohfields(i)%ref%boundary_comm_infos
                do j = 1, size(infos)
                    id = infos(j)%id
                    self%boundary_communication_infos(id) = infos(j)
                end do
            end do
        end block
    end subroutine

    subroutine ohhelp_allocate_ohfield(self, ohfield)
        class(t_OhHelp), intent(in) :: self
        class(t_OhField), intent(inout) :: ohfield

        integer :: nelements
        integer :: nfields
        integer :: field_size(2, 3)

        nelements = ohfield%extension_info%nelements
        nfields = ohfield%nfields

        field_size(:, :) = self%field_sizes(:, :, ohfield%extension_info%id)

        allocate (ohfield%values(nelements, &
                                 field_size(1, 1):field_size(2, 1), &
                                 field_size(1, 2):field_size(2, 2), &
                                 field_size(1, 3):field_size(2, 3), &
                                 nfields))
    end subroutine

    subroutine ohhelp_correct_load_balancing(self, ohparticles, eb, ohfields_to_be_notified)
        class(t_OhHelp), intent(inout) :: self
        class(t_OhParticles), intent(inout) :: ohparticles
        type(tp_OhField), intent(inout) :: ohfields_to_be_notified(:)
        class(t_OhField), intent(inout) :: eb

        call self%transbound(ohparticles)

        if (.not. self%requires_broadcast_field()) then
            return
        end if

        call self%broadcast_field(eb)
        self%current_mode = OHH_SECONDARY_MODE

        block
            integer :: i

            do i = 1, size(ohfields_to_be_notified)
                call self%notify_subdomain_id(ohfields_to_be_notified(i)%ref)
            end do
        end block
    end subroutine

    subroutine ohhelp_transbound(self, ohparticles)
        class(t_OhHelp), intent(inout) :: self
        class(t_OhParticles), intent(inout) :: ohparticles
        integer :: status

        self%current_mode = oh3_transbound(self%current_mode, status)

        ohparticles%particle_count_histgram(self%subdomain_id(1) + 1, :, 1) &
            = ohparticles%total_local_particles(:, 1)

        if (self%subdomain_id(2) >= 0) then
            ohparticles%particle_count_histgram(self%subdomain_id(2) + 1, :, 2) &
                = ohparticles%total_local_particles(:, 2)
        end if
    end subroutine

    function ohhelp_requires_broadcast_field(self) result(ret)
        class(t_OhHelp), intent(inout) :: self
        logical :: ret

        ret = self%current_mode == OHH_REQUIRES_BCAST
    end function

    subroutine ohhelp_broadcast_field(self, ohfield)
        class(t_OhHelp), intent(inout) :: self
        type(t_OhField), intent(inout) :: ohfield

        call oh3_bcast_field(ohfield%values(1, 0, 0, 0, 1), &
                             ohfield%values(1, 0, 0, 0, 2), &
                             ohfield%extension_info%id)
    end subroutine

    subroutine check_particles_in_subdomain(particles, pbase, primary_or_secondary)
        type(oh_particle), intent(inout) :: particles(:)
        integer, intent(in) :: pbase(3)
        integer, intent(in) :: primary_or_secondary
        type(oh_particle) :: particle

        integer :: nid

        nid = oh3_map_particle_to_neighbor(particle%x, particle%y, particle%z, primary_or_secondary)
    end subroutine

    subroutine ohhelp_notify_subdomain_id(self, ohfield)
        class(t_OhHelp), intent(inout) :: self
        type(t_OhField), intent(inout) :: ohfield

        ohfield%subdomain_range(:, :, 1) = self%subdomain_range(:, :, self%subdomain_id(1)+1)

        if (self%is_secondary_mode()) then
            ohfield%subdomain_range(:, :, 2) = self%subdomain_range(:, :, self%subdomain_id(2)+1)
        end if
    end subroutine

    subroutine ohhelp_inject_particle(self, particle)
        class(t_OhHelp), intent(inout) :: self
        type(oh_particle), intent(in) :: particle

        call oh2_inject_particle(particle)
    end subroutine

    function ohhelp_is_primary_mode(self) result(ret)
        class(t_OhHelp), intent(inout) :: self
        logical :: ret

        ret = self%current_mode == OHH_PRIMARY_MODE
    end function

    function ohhelp_is_secondary_mode(self) result(ret)
        class(t_OhHelp), intent(inout) :: self
        logical :: ret

        ret = self%current_mode == OHH_SECONDARY_MODE
    end function

end module
