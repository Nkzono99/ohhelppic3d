
module m_block_communicator_list
    use m_list, only: t_List, init_list
    use m_block_communicator, only: t_BlockCommunicator
    implicit none

    private
    public t_BlockCommunicatorList
    public new_BlockCommunicatorList

    type, extends(t_List) :: t_BlockCommunicatorList
        type(t_BlockCommunicator), allocatable :: buffer(:)
        type(t_BlockCommunicator), allocatable :: tmp_buffer(:)

    contains
        procedure :: allocate_tmp_buffer => blockSendReceiverList_allocate_tmp_buffer
        procedure :: copy_to_tmp_buffer => blockSendReceiverList_copy_to_tmp_buffer
        procedure :: switch_to_tmp_buffer => blockSendReceiverList_switch_to_tmp_buffer
        procedure :: destroy => blockSendReceiverList_destroy
        procedure :: append => blockSendReceiverList_append
        procedure :: get => blockSendReceiverList_get
    end type

contains

    function new_BlockCommunicatorList(max_size, growth_factor) result(obj)
        integer, intent(in), optional :: max_size
        double precision, intent(in), optional :: growth_factor
        type(t_BlockCommunicatorList) :: obj

        call init_list(obj, max_size, growth_factor)
    end function

    subroutine blockSendReceiverList_allocate_tmp_buffer(self, n)
        class(t_BlockCommunicatorList), intent(inout) :: self
        integer, intent(in) :: n

        allocate (self%tmp_buffer(n))
    end subroutine

    subroutine blockSendReceiverList_copy_to_tmp_buffer(self)
        class(t_BlockCommunicatorList), intent(inout) :: self

        integer :: i

        do i = 1, self%current_size
            self%tmp_buffer(i) = self%buffer(i)
        end do

        deallocate (self%buffer)
    end subroutine

    subroutine blockSendReceiverList_switch_to_tmp_buffer(self)
        class(t_BlockCommunicatorList), intent(inout) :: self

        allocate (self%buffer, source=self%tmp_buffer)
        deallocate (self%tmp_buffer)
    end subroutine

    subroutine blockSendReceiverList_append(self, setting)
        class(t_BlockCommunicatorList), intent(inout) :: self
        type(t_BlockCommunicator), intent(in) :: setting

        if (self%current_size == self%max_size) then
            call self%extent_size()
        end if

        self%current_size = self%current_size + 1
        self%buffer(self%current_size) = setting
    end subroutine

    subroutine blockSendReceiverList_destroy(self)
        class(t_BlockCommunicatorList), intent(inout) :: self

        type(t_BlockCommunicator) :: setting
        integer :: i

        do i = 1, self%current_size
            setting = self%get(i)
            call setting%destroy
        end do

        deallocate (self%buffer)
    end subroutine

    function blockSendReceiverList_get(self, i) result(ret)
        class(t_BlockCommunicatorList), intent(in) :: self
        integer, intent(in) :: i
        type(t_BlockCommunicator) :: ret

        ret = self%buffer(i)
    end function

end module
