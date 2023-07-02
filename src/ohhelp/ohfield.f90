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
        integer :: id !> Identifier
        integer :: nelements = 0 !> The number of elements.
        integer :: nextensions(2) !>  The size of subdomain expansion (lower, upper)
        integer :: nextensions_for_broadcast(2) ! The size of subdomain expansion when broadcast (lower, upper)
        integer :: nextensions_for_reduction(2) ! The size of subdomain expansion when reduction (lower, upper)
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
        integer :: id !> boundary communication type identifier
        type(t_BoundaryCommunicationInfo) :: infos(NBOUNDARY_CONDITION_TYPES)
    end type

    type t_OhField
        double precision, allocatable :: values(:, :, :, :, :) ! (nelements, nlx, nly, nlz, nfields)

        type(t_FieldExtensionInfo) :: extension_info

        integer :: nfields
        integer :: field_size(2, 3)

        integer :: nboundary_comm_infos
        type(t_BoundaryCommunicationInfos), allocatable :: boundary_comm_infos(:)
    contains
        procedure :: make_copy => ohfield_make_copy
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
        integer, intent(in) :: id !> boundary communication type identifier
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
        obj%nfields = nfields

        if (present(boundary_comm_infos)) then
            obj%nboundary_comm_infos = size(boundary_comm_infos)

            allocate (obj%boundary_comm_infos(obj%nboundary_comm_infos))
            obj%boundary_comm_infos(:) = boundary_comm_infos(:)
        else
            obj%nboundary_comm_infos = 0
        end if
    end function

    function ohfield_make_copy(self, nfields) result(copy)
        class(t_OhField), intent(in) :: self
        integer, intent(in) :: nfields

        type(t_OhField) :: copy

        copy = new_OhField(self%extension_info, nfields, self%boundary_comm_infos)
    end function

end module
