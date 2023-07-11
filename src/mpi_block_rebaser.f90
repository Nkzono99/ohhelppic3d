! TODO: モジュール名とクラス名をよりわかりやすい名前に変える(rebaseはよくわからない)
module m_mpi_block_rebaser
    use m_block
    use m_block_list
    use m_block_communicator
    use m_block_communicator_list
    use mpi
    use m_get_default, only: get_default
    implicit none

    private
    public t_MPIBlockRebaser
    public new_MPIBlockRebaser

    type t_MPIBlockRebaser
        !> Block that local process has.
        type(t_Block) :: local_block
        !> Block that local process requires.
        type(t_Block) :: require_block

        type(t_BlockCommunicatorList), private :: block_senders
        type(t_BlockCommunicatorList), private :: block_receivers

        integer, allocatable :: pids(:)
        integer :: ipid

        integer(kind=kind(MPI_COMM_WORLD)), private :: comm
    contains
        procedure :: rebase => mpiBlockRebaser_rebase
        procedure :: destroy => mpiBlockRebaser_destroy
    end type

    interface new_MPIBlockRebaser
        procedure :: new_MPIBlockRebaser_with_blocks
        procedure :: new_MPIBlockRebaser_with_local_block
    end interface

contains

    function new_MPIBlockRebaser_with_blocks(local_blocks, require_blocks, pids, ipid, comm) result(obj)
        !> Blocks that each process has.
        type(t_BlockList), intent(in) :: local_blocks
        !> Blocks that each process requires.
        type(t_BlockList), intent(in) :: require_blocks
        !> Array of process ID (= rank = 0, 1, 2, ..., nproc-1).
        !> It should be the same as the order of blocks and require_blocks.
        integer, intent(in) :: pids(:)
        !> Index to identify the current process ID from pids(1:len(pids)).
        integer, intent(in) :: ipid
        !> MPI Communicator.
        integer(kind=kind(MPI_COMM_WORLD)), intent(in) :: comm
        type(t_MPIBlockRebaser) :: obj

        obj%local_block = local_blocks%get(ipid)
        obj%require_block = require_blocks%get(ipid)

        obj%pids = pids
        obj%ipid = ipid

        obj%block_senders = new_BlockCommunicatorList()
        obj%block_receivers = new_BlockCommunicatorList()

        block ! Send settings
            type(t_Block) :: overlapped
            type(t_BlockCommunicator) :: sender

            integer :: i

            do i = 1, require_blocks%current_size
                overlapped = obj%local_block%overlapped(require_blocks%get(i))

                if (overlapped%size() == 0) then
                    cycle
                end if

                sender = new_BlockCommunicator(obj%local_block, overlapped, pids(i), comm)
                call obj%block_senders%append(sender)
            end do
        end block

        block ! Recv settings
            type(t_Block) :: overlapped
            type(t_BlockCommunicator) :: receiver

            integer :: i

            do i = 1, local_blocks%current_size
                overlapped = obj%require_block%overlapped(local_blocks%get(i))

                if (overlapped%size() == 0) then
                    cycle
                end if

                receiver = new_BlockCommunicator(obj%require_block, overlapped, pids(i), comm)
                call obj%block_receivers%append(receiver)
            end do
        end block

        obj%comm = comm
    end function

    function new_MPIBlockRebaser_with_local_block(local_block, require_block, pids, ipid, comm, tag) result(obj)
        !> Blocks that each process has.
        type(t_Block), intent(in) :: local_block
        !> Blocks that each process requires.
        type(t_Block), intent(in) :: require_block
        !> Array of process ID (= rank = 0, 1, 2, ..., nproc-1).
        !> It should be the same as the order of blocks and require_blocks.
        integer, intent(in) :: pids(:)
        !> Index to identify the current process ID from pids(1:len(pids)).
        integer, intent(in) :: ipid
        !> MPI Communicator.
        integer(kind=kind(MPI_COMM_WORLD)), intent(in) :: comm
        integer, optional, intent(in) :: tag
        type(t_MPIBlockRebaser) :: obj

        !> Blocks that each process has.
        type(t_BlockList) :: local_blocks
        !> Blocks that each process requires.
        type(t_BlockList) :: require_blocks

        local_blocks = mpi_collect_blocks(local_block, pids, comm, tag)
        require_blocks = mpi_collect_blocks(require_block, pids, comm, tag)

        obj = new_MPIBlockRebaser(local_blocks, require_blocks, pids, ipid, comm)
    end function

    function mpi_collect_blocks(block, pids, comm, tag) result(blocks)
        type(t_Block), intent(in) :: block
        !> Array of process ID (= rank = 0, 1, 2, ..., nproc-1).
        !> It should be the same as the order of blocks and require_blocks.
        integer, intent(in) :: pids(:)
        integer, intent(in) :: comm
        integer, optional, intent(in) :: tag
        type(t_BlockList) :: blocks

        integer :: send_data(SIZE_OF_BLOCK_ARRAY)
        integer :: recv_datas(SIZE_OF_BLOCK_ARRAY, size(pids))
        integer :: send_requests(size(pids))
        integer :: recv_requests(size(pids))
        integer :: send_status(MPI_STATUS_SIZE, size(pids))
        integer :: recv_status(MPI_STATUS_SIZE, size(pids))
        integer :: ierr

        integer :: i, ip

        call block%to_array(send_data(:))

        do i = 1, size(pids)
            ip = pids(i)
            call MPI_Isend(send_data(1), SIZE_OF_BLOCK_ARRAY, MPI_INTEGER, &
                           ip, &
                           get_default(tag, 0), comm, &
                           send_requests(i), ierr)
            call MPI_Irecv(recv_datas(1, i), SIZE_OF_BLOCK_ARRAY, MPI_INTEGER, &
                           ip, &
                           get_default(tag, 0), comm, &
                           recv_requests(i), ierr)
        end do

        call MPI_Waitall(size(pids), recv_requests(:), recv_status(:, :), ierr)
        call MPI_Waitall(size(pids), send_requests(:), send_status(:, :), ierr)

        blocks = new_BlockList()
        block
            type(t_Block) :: blk
            do i = 1, size(pids)
                call blk%from_array(recv_datas(:, i))
                call blocks%append(blk)
            end do
        end block
    end function

    subroutine mpiBlockRebaser_rebase(self, send_data, recv_data, tag)
        class(t_MPIBlockRebaser), intent(in) :: self
        double precision, intent(in) :: send_data(self%local_block%start(1):self%local_block%end(1), &
                                                  self%local_block%start(2):self%local_block%end(2), &
                                                  self%local_block%start(3):self%local_block%end(3))
        double precision, intent(inout) :: recv_data(self%require_block%start(1):self%require_block%end(1), &
                                                     self%require_block%start(2):self%require_block%end(2), &
                                                     self%require_block%start(3):self%require_block%end(3))
        integer, intent(in), optional :: tag

        integer :: send_requests(self%block_senders%current_size)
        integer :: recv_requests(self%block_receivers%current_size)
        integer :: send_status(MPI_STATUS_SIZE, self%block_senders%current_size)
        integer :: recv_status(MPI_STATUS_SIZE, self%block_receivers%current_size)

        block
            integer :: i
            type(t_BlockCommunicator) :: sender

            do i = 1, self%block_senders%current_size
                sender = self%block_senders%get(i)
                ! print *, self%pids(self%ipid), 'send', sender%pid
                call sender%isend(send_data(:, :, :), get_default(tag, 0), send_requests(i))
            end do
        end block

        block
            integer :: i
            type(t_BlockCommunicator) :: receiver
            do i = 1, self%block_receivers%current_size
                receiver = self%block_receivers%get(i)
                ! print *, self%pids(self%ipid), 'recv', receiver%pid
                call receiver%irecv(recv_data(:, :, :), get_default(tag, 0), recv_requests(i))
            end do
        end block

        ! print *, self%pids(self%ipid), self%block_senders%current_size, self%block_receivers%current_size

        block
            integer :: ierr
            call MPI_Waitall(self%block_receivers%current_size, recv_requests(:), recv_status(:, :), ierr)

            call MPI_Waitall(self%block_senders%current_size, send_requests(:), send_status(:, :), ierr)
        end block
    end subroutine

    subroutine mpiBlockRebaser_destroy(self)
        class(t_MPIBlockRebaser), intent(inout) :: self

        call self%block_senders%destroy
        call self%block_receivers%destroy
    end subroutine

end module
