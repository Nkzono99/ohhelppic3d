module m_ohfield_factory
    use m_ohfield
    implicit none

    private
    public t_OhFieldFactory
    public new_OhFieldFactory

    type t_OhFieldFactory
        integer :: extension_id_count
        integer :: comm_id_count
    contains
        procedure :: create_field => factory_create_field
        procedure, private :: create_electromagnetic_field => factory_create_electromagnetic_field
        procedure, private :: create_current_field => factory_create_current_field
        procedure, private :: create_density_field => factory_create_density_field
        procedure, private :: create_potential_field => factory_create_potential_field
        procedure, private :: create_field_wrapper => factory_create_field_wrapper
    end type

contains

    function new_OhFieldFactory() result(obj)
        type(t_OhFieldFactory) :: obj

        obj%extension_id_count = 1
        obj%comm_id_count = 1
    end function

    function factory_create_field(self, name, nfields) result(field)
        class(t_OhFieldFactory), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer, optional, intent(in) :: nfields
        type(t_OhField) :: field

        integer :: nfields_ = 2

        if (present(nfields)) then
            nfields_ = nfields
        end if

        select case (name)
        case ('eb', 'electromagnetic')
            field = self%create_electromagnetic_field(nfields_)
        case ('j', 'current')
            field = self%create_current_field(nfields_)
        case ('rho', 'density')
            field = self%create_density_field(nfields_)
        case ('phi', 'potential')
            field = self%create_potential_field(nfields_)
        case default
            print *, 'InvalidArgumentError: input name is invalid (OhFieldFactory%create_field):', name
            return
        end select
    end function

    function factory_create_field_wrapper(self, &
                                          nelements, &
                                          nfields, &
                                          nextensions, &
                                          nextensions_for_broadcast, &
                                          nextensions_for_reduction, &
                                          downward_comm_send_offset, &
                                          downward_comm_receive_offset, &
                                          downward_comm_nsends, &
                                          upward_comm_send_offset, &
                                          upward_comm_receive_offset, &
                                          upward_comm_nsends) result(ohfield)
        class(t_OhFieldFactory), intent(inout) :: self
        integer, intent(in) :: nelements
        integer, intent(in) :: nfields
        integer, intent(in) :: nextensions(2)
        integer, intent(in) :: nextensions_for_broadcast(2)
        integer, intent(in) :: nextensions_for_reduction(2)
        integer, intent(in) :: downward_comm_send_offset
        integer, intent(in) :: downward_comm_receive_offset
        integer, intent(in) :: downward_comm_nsends
        integer, intent(in) :: upward_comm_send_offset
        integer, intent(in) :: upward_comm_receive_offset
        integer, intent(in) :: upward_comm_nsends

        type(t_OhField) :: ohfield

        type(t_FieldExtensionInfo) :: exinfo
        type(t_BoundaryCommunicationInfos) :: bcomminfos(1)

        exinfo = new_FieldExtensionInfo(id=self%extension_id_count, &
                                        nelements=nelements, &
                                        nextensions=nextensions, &
                                        nextensions_for_broadcast=nextensions_for_broadcast, &
                                        nextensions_for_reduction=nextensions_for_reduction)
        self%extension_id_count = self%extension_id_count + 1

        block
            type(t_BoundaryCommunicationInfo) :: infos(NBOUNDARY_CONDITION_TYPES)

            infos(1) = new_BoundaryCommunicationInfo( &
                       downward_comm_send_offset, &
                       downward_comm_receive_offset, &
                       downward_comm_nsends, &
                       upward_comm_send_offset, &
                       upward_comm_receive_offset, &
                       upward_comm_nsends)

            infos(2) = new_BoundaryCommunicationInfo( &
                       downward_comm_send_offset, &
                       downward_comm_receive_offset, &
                       downward_comm_nsends, &
                       upward_comm_send_offset, &
                       upward_comm_receive_offset, &
                       upward_comm_nsends)

            bcomminfos(1) = &
                new_BoundaryCommunicationInfos(self%comm_id_count, infos)

            self%comm_id_count = self%comm_id_count + 1
        end block

        ohfield = new_OhField(exinfo, nfields, bcomminfos)
    end function

    function factory_create_electromagnetic_field(self, nfields) result(field)
        class(t_OhFieldFactory), intent(inout) :: self
        integer, intent(in) :: nfields
        type(t_OhField) :: field

        field = self%create_field_wrapper( &
                nelements=6, &
                nfields=nfields, &
                nextensions=[-1, 2], &
                nextensions_for_broadcast=[-1, 2], &
                nextensions_for_reduction=[0, 0], &
                downward_comm_send_offset=0, &
                downward_comm_receive_offset=0, &
                downward_comm_nsends=2, &
                upward_comm_send_offset=-1, &
                upward_comm_receive_offset=-1, &
                upward_comm_nsends=1)
    end function

    function factory_create_current_field(self, nfields) result(field)
        class(t_OhFieldFactory), intent(inout) :: self
        integer, intent(in) :: nfields
        type(t_OhField) :: field

        field = self%create_field_wrapper( &
                nelements=3, &
                nfields=nfields, &
                nextensions=[0, 0], &
                nextensions_for_broadcast=[0, 0], &
                nextensions_for_reduction=[-1, 2], &
                downward_comm_send_offset=-1, &
                downward_comm_receive_offset=2, &
                downward_comm_nsends=3, &
                upward_comm_send_offset=-1, &
                upward_comm_receive_offset=-4, &
                upward_comm_nsends=3)
    end function

    function factory_create_density_field(self, nfields) result(field)
        class(t_OhFieldFactory), intent(inout) :: self
        integer, intent(in) :: nfields
        type(t_OhField) :: field

        field = self%create_field_wrapper( &
                nelements=1, &
                nfields=nfields, &
                nextensions=[0, 1], &
                nextensions_for_broadcast=[0, 0], &
                nextensions_for_reduction=[0, 1], &
                downward_comm_send_offset=0, &
                downward_comm_receive_offset=1, &
                downward_comm_nsends=1, &
                upward_comm_send_offset=0, &
                upward_comm_receive_offset=-1, &
                upward_comm_nsends=1)
    end function

    function factory_create_potential_field(self, nfields) result(field)
        class(t_OhFieldFactory), intent(inout) :: self
        integer, intent(in) :: nfields
        type(t_OhField) :: field

        field = self%create_field_wrapper( &
                nelements=1, &
                nfields=nfields, &
                nextensions=[0, 0], &
                nextensions_for_broadcast=[-1, 2], &
                nextensions_for_reduction=[0, 0], &
                downward_comm_send_offset=0, &
                downward_comm_receive_offset=0, &
                downward_comm_nsends=1, &
                upward_comm_send_offset=-1, &
                upward_comm_receive_offset=-1, &
                upward_comm_nsends=1)
    end function

end module
