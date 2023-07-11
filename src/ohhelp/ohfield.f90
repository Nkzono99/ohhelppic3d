module m_ohfield
    implicit none

    private
    public NBOUNDARY_CONDITION_TYPES
    public BOUNDARY_CONDITION_PERIODIC
    public BOUNDARY_CONDITION_NO_PERIODIC
    public t_BoundaryCommunicationInfo, new_BoundaryCommunicationInfo
    public t_BoundaryCommunicationInfos, new_BoundaryCommunicationInfos
    public t_FieldExtensionInfo, new_FieldExtensionInfo
    public t_OhField, new_OhField
    public tp_OhField

    integer, parameter :: NBOUNDARY_CONDITION_TYPES = 2
    integer, parameter :: BOUNDARY_CONDITION_PERIODIC = 1
    integer, parameter :: BOUNDARY_CONDITION_NO_PERIODIC = 2

    type t_FieldExtensionInfo
        !> Identifier
        integer :: id
        !> The number of elements.
        integer :: nelements = 0
        !>  The size of subdomain expansion [lower, upper]
        integer :: nextensions(2)
        !> The size of subdomain expansion when broadcast [lower, upper]
        integer :: nextensions_for_broadcast(2)
        !> The size of subdomain expansion when reduction [lower, upper]
        integer :: nextensions_for_reduction(2)
    end type

    type t_BoundaryCommunicationInfo
        !> Offset to send at downward communication
        integer :: downward_comm_send_offset
        !> Offset to receive at downward communication
        integer :: downward_comm_receive_offset
        !> The number of planes to send/receive at downward communication
        integer :: downward_comm_nsends

        !> Offset to send at upward communication
        integer :: upward_comm_send_offset
        !> Offset to receive at upward communication
        integer :: upward_comm_receive_offset
        !> The number of planes to send/receive at upward communication
        integer :: upward_comm_nsends
    end type

    type t_BoundaryCommunicationInfos
        !> boundary communication type identifier
        integer :: id
        type(t_BoundaryCommunicationInfo) :: infos(NBOUNDARY_CONDITION_TYPES)
    end type

    type t_OhField
        double precision, allocatable :: values(:, :, :, :, :) ! (nelements, nlx, nly, nlz, nfields)

        type(t_FieldExtensionInfo) :: extension_info

        integer :: nelements
        integer :: nfields
        integer :: field_size(2, 3)
        integer :: subdomain_range(2, 3, 2)

        integer :: nboundary_comm_infos
        type(t_BoundaryCommunicationInfos), allocatable :: boundary_comm_infos(:)
    contains
        procedure :: make_copy => ohfield_make_copy
        procedure :: allocate => ohfield_allocate
        procedure :: to_global_index => ohfield_to_global_index
        procedure :: to_local_index => ohfield_to_local_index
        procedure :: to_local_position => ohfield_to_local_position
        procedure :: local_shape => ohfield_local_shape
    end type

    type tp_OhField
        class(t_OhField), pointer :: ref
    end type

contains

    function new_FieldExtensionInfo(id, &
                                    nelements, &
                                    nextensions, &
                                    nextensions_for_broadcast, &
                                    nextensions_for_reduction) result(obj)
        integer, intent(in) :: id
        integer, intent(in) :: nelements
        integer, intent(in) :: nextensions(2)
        integer, intent(in) :: nextensions_for_broadcast(2)
        integer, intent(in) :: nextensions_for_reduction(2)

        type(t_FieldExtensionInfo) :: obj

        obj%id = id
        obj%nelements = nelements
        obj%nextensions = nextensions
        obj%nextensions_for_broadcast = nextensions_for_broadcast
        obj%nextensions_for_reduction = nextensions_for_reduction
    end function

    function new_BoundaryCommunicationInfo(downward_comm_send_offset, &
                                           downward_comm_receive_offset, &
                                           downward_comm_nsends, &
                                           upward_comm_send_offset, &
                                           upward_comm_receive_offset, &
                                           upward_comm_nsends) result(obj)
        !> Offset to send at downward communication
        integer, intent(in) :: downward_comm_send_offset
        !> Offset to receive at downward communication
        integer, intent(in) :: downward_comm_receive_offset
        !> The number of planes to send/receive at downward communication
        integer, intent(in) :: downward_comm_nsends
        !> Offset to send at upward communication
        integer, intent(in) :: upward_comm_send_offset
        !> Offset to receive at upward communication
        integer, intent(in) :: upward_comm_receive_offset
        !> The number of planes to send/receive at upward communication
        integer, intent(in) :: upward_comm_nsends

        type(t_BoundaryCommunicationInfo) :: obj

        obj%downward_comm_send_offset = downward_comm_send_offset
        obj%downward_comm_receive_offset = downward_comm_receive_offset
        obj%downward_comm_nsends = downward_comm_nsends
        obj%upward_comm_send_offset = upward_comm_send_offset
        obj%upward_comm_receive_offset = upward_comm_receive_offset
        obj%upward_comm_nsends = upward_comm_nsends
    end function

    function new_BoundaryCommunicationInfos(id, boundary_communication_infos) result(obj)
        !> boundary communication type identifier
        integer, intent(in) :: id
        type(t_BoundaryCommunicationInfo) :: boundary_communication_infos(NBOUNDARY_CONDITION_TYPES)
        type(t_BoundaryCommunicationInfos) :: obj

        obj%id = id
        obj%infos(:) = boundary_communication_infos(:)
    end function

    function new_OhField(extension_info, nfields, boundary_comm_infos) result(obj)
        type(t_FieldExtensionInfo), intent(in) :: extension_info
        integer, intent(in) :: nfields
        type(t_BoundaryCommunicationInfos), intent(in), optional :: boundary_comm_infos(:)
        type(t_OhField) :: obj

        obj%extension_info = extension_info
        obj%nelements = extension_info%nelements
        obj%nfields = nfields

        if (present(boundary_comm_infos)) then
            obj%nboundary_comm_infos = size(boundary_comm_infos)

            allocate (obj%boundary_comm_infos(obj%nboundary_comm_infos))
            obj%boundary_comm_infos(:) = boundary_comm_infos(:)
        else
            obj%nboundary_comm_infos = 0
        end if
    end function

    subroutine ohfield_allocate(self, nelements, field_size, nfields)
        class(t_Ohfield), intent(inout) :: self
        integer, intent(in) :: nelements
        integer, intent(in) :: field_size(2, 3)
        integer, intent(in) :: nfields

        allocate (self%values(nelements, &
                              field_size(1, 1):field_size(2, 1), &
                              field_size(1, 2):field_size(2, 2), &
                              field_size(1, 3):field_size(2, 3), &
                              nfields))
    end subroutine

    function ohfield_make_copy(self, nfields) result(copy)
        class(t_OhField), intent(in) :: self
        integer, intent(in) :: nfields

        type(t_OhField) :: copy

        copy = new_OhField(self%extension_info, nfields, self%boundary_comm_infos)
    end function

    function ohfield_to_global_index(self, local_index, ps) result(global_index)
        class(t_Ohfield), intent(in) :: self
        integer, intent(in) :: local_index(3)
        integer, intent(in) :: ps
        integer :: global_index(3)

        global_index = self%subdomain_range(1, :, ps) + local_index(:)
    end function

    function ohfield_to_local_index(self, global_index, ps) result(local_index)
        class(t_Ohfield), intent(in) :: self
        integer, intent(in) :: global_index(3)
        integer, intent(in) :: ps
        integer :: local_index(3)

        local_index = global_index - self%subdomain_range(1, :, ps)
    end function

    function ohfield_to_local_position(self, global_position, ps) result(ret)
        class(t_Ohfield), intent(in) :: self
        double precision, intent(in) :: global_position(3)
        integer, intent(in) :: ps
        double precision :: ret(3)

        ret = global_position - self%subdomain_range(1, :, ps)
    end function

    function ohfield_local_shape(self, ps) result(ret)
        class(t_Ohfield), intent(in) :: self
        integer, intent(in) :: ps
        integer :: ret(3)

        ret = self%subdomain_range(2, :, ps) - self%subdomain_range(1, :, ps) + 1
    end function
end module
