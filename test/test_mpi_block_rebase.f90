module test_mpi_block_rebaser
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use m_mpi_block_rebaser
    use m_block
    use m_check_array, only: check
    use mpi
    use m_str
    use m_mpi_block
    implicit none

    private
    public collect_mpi_block_rebaser

    integer :: nprocs, myid

contains

    subroutine collect_mpi_block_rebaser(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        integer :: ierr

        call mpi_init(ierr)
        if (ierr /= 0) then
            write (error_unit, '(a)') 'Use "mpiexec -n 4" to run this test.'
            stop 1
        end if

        call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
        call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)

        testsuite = [ &
                    new_unittest("rebase", test_rebase) &
                    ]
    end subroutine

    subroutine test_rebase(error)
        type(error_type), allocatable, intent(out) :: error

        type(t_Block) :: global_block
        type(t_Block) :: local_block
        type(t_Block) :: require_block

        integer :: n(3)
        integer :: nnodes(3)
        integer :: ids(3)

        n = [16, 32, 64]
        nnodes = [2, 1, 2]
        global_block = new_Block([1, 1, 1], n)
        call create_local_block(global_block, nnodes, myid, local_block, ids)
        require_block = new_Block([1, 1, ids*n(3)/4 + 1], [n(1), n(2), (ids + 1)*n(3)/4])

        call check_rebase(local_block, require_block)
        if (allocated(error)) return

        n = [16, 32, 64]
        nnodes = [1, 2, 2]
        global_block = new_Block([1, 1, 1], n)
        call create_local_block(global_block, nnodes, myid, local_block, ids)
        require_block = new_Block([1, 1, ids*n(3)/4 + 1], [n(1), n(2), (ids + 1)*n(3)/4])

        call check_rebase(local_block, require_block)
        if (allocated(error)) return

        n = [16, 32, 64]
        nnodes = [4, 1, 1]
        global_block = new_Block([1, 1, 1], n)
        call create_local_block(global_block, nnodes, myid, local_block, ids)
        require_block = new_Block([1, 1, ids*n(3)/4 + 1], [n(1), n(2), (ids + 1)*n(3)/4])

        call check_rebase(local_block, require_block)
        if (allocated(error)) return

    contains

        subroutine check_rebase(lblk, rblk)
            type(t_Block), intent(in) :: lblk, rblk

            type(t_MpiBlockRebaser) :: rebaser

            double precision :: send_data(lblk%start(1):lblk%end(1), &
                                          lblk%start(2):lblk%end(2), &
                                          lblk%start(3):lblk%end(3))
            double precision ::  recv_data(rblk%start(1):rblk%end(1), &
                                           rblk%start(2):rblk%end(2), &
                                           rblk%start(3):rblk%end(3))
            integer :: i, j, k

            double precision :: expected(rblk%start(1):rblk%end(1), &
                                         rblk%start(2):rblk%end(2), &
                                         rblk%start(3):rblk%end(3))

            do concurrent(i=lblk%start(1):lblk%end(1), j=lblk%start(2):lblk%end(2), k=lblk%start(3):lblk%end(3))
                send_data(i, j, k) = i + 2*j + 3*k
            end do
            do concurrent(i=rblk%start(1):rblk%end(1), j=rblk%start(2):rblk%end(2), k=rblk%start(3):rblk%end(3))
                expected(i, j, k) = i + 2*j + 3*k
            end do

            rebaser = new_MPIBlockRebaser(lblk, rblk, [0, 1, 2, 3], myid + 1, MPI_COMM_WORLD)

            call rebaser%rebase(send_data, recv_data)

            call check(error, recv_data, expected, th=1d-3)

            call rebaser%destroy()
        end subroutine
    end subroutine
end module
