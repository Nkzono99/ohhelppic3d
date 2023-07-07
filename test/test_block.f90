module test_block
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use m_block
    use m_check_array, only: check
    implicit none

    private
    public collect_block

contains

    subroutine collect_block(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("new_Block", test_new_Block), &
                    new_unittest("block.size", test_size), &
                    new_unittest("block.to_array", test_to_array), &
                    new_unittest("block.from_array", test_from_array), &
                    new_unittest("block.overlapped", test_overlapped) &
                    ]
    end subroutine

    subroutine test_new_Block(error)
        type(error_type), allocatable, intent(out) :: error

        call check_start_end([0, 0, 0], [0, 0, 0])
        if (allocated(error)) return

        call check_start_end([1, 2, 3], [4, 5, 6])
        if (allocated(error)) return

        call check_start_end([-2, -3, -5], [4, 5, 6])
        if (allocated(error)) return

        call check_start_end([4, 5, 6], [-2, -3, -5])
        if (allocated(error)) return

    contains

        subroutine check_start_end(s, e)
            integer, intent(in) :: s(3), e(3)
            type(t_Block) :: blk

            integer :: i

            blk = new_Block(s, e)

            call check(error, blk%start, s)
            call check(error, blk%end, e)
        end subroutine
    end subroutine

    subroutine test_size(error)
        type(error_type), allocatable, intent(out) :: error

        ! The block defines a closed interval ([start, end]), so its size is 1.
        call check_size([0, 0, 0], [0, 0, 0], 1)
        if (allocated(error)) return

        call check_size([0, 0, 0], [1, 0, 0], 2)
        if (allocated(error)) return

        call check_size([0, 0, 0], [0, 0, 1], 2)
        if (allocated(error)) return

        call check_size([0, 0, 0], [0, 1, 0], 2)
        if (allocated(error)) return

        call check_size([1, 2, 3], [4, 7, 10], 4*6*8)
        if (allocated(error)) return

        ! Size of the block is 0 when start > end.
        call check_size([4, 7, 10], [1, 2, 3], 0)
        if (allocated(error)) return

        call check_size([0, 0, 0], [-1, 0, 0], 0)
        if (allocated(error)) return

        call check_size([0, 0, 0], [0, -1, 0], 0)
        if (allocated(error)) return

        call check_size([0, 0, 0], [0, 0, -1], 0)
        if (allocated(error)) return

    contains
        subroutine check_size(s, e, sz)
            integer, intent(in) :: s(3), e(3), sz
            type(t_Block) :: blk

            blk = new_Block(s, e)
            call check(error, blk%size(), sz)
            if (allocated(error)) return
        end subroutine
    end subroutine

    subroutine test_overlapped(error)
        type(error_type), allocatable, intent(out) :: error

        call check_overlapped([0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0])
        if (allocated(error)) return

        call check_overlapped([0, 0, 0], [10, 11, 12], [3, 4, 5], [6, 7, 8], [3, 4, 5], [6, 7, 8])
        if (allocated(error)) return

        call check_overlapped([4, 5, 6], [10, 11, 12], [3, 4, 5], [6, 7, 8], [4, 5, 6], [6, 7, 8])
        if (allocated(error)) return

        call check_overlapped([0, 0, 0], [4, 5, 6], [3, 4, 5], [6, 7, 8], [3, 4, 5], [4, 5, 6])
        if (allocated(error)) return

        call check_overlapped([0, 0, 0], [0, -1, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [-1, -1, -1])
        if (allocated(error)) return

        call check_overlapped([0, 0, 0], [0, 0, 0], [0, 0, 0], [0, -1, 0], [0, 0, 0], [-1, -1, -1])
        if (allocated(error)) return

        call check_overlapped([0, 0, 0], [0, -1, 0], [0, 0, 0], [0, -1, 0], [0, 0, 0], [-1, -1, -1])
        if (allocated(error)) return

    contains

        subroutine check_overlapped(s1, e1, s2, e2, s3, e3)
            integer, intent(in) :: s1(3), e1(3)
            integer, intent(in) :: s2(3), e2(3)
            integer, intent(in) :: s3(3), e3(3)
            type(t_Block) :: blk1, blk2, overlapped

            blk1 = new_Block(s1, e1)
            blk2 = new_Block(s2, e2)

            overlapped = blk1%overlapped(blk2)

            call check(error, overlapped%start, s3)
            if (allocated(error)) return

            call check(error, overlapped%end, e3)
            if (allocated(error)) return
        end subroutine
    end subroutine

    subroutine test_to_array(error)
        type(error_type), allocatable, intent(out) :: error

        call check_to_array([0, 0, 0], [0, 0, 0], [0, 0, 0, 0, 0, 0])
        if (allocated(error)) return

        call check_to_array([1, 2, 3], [4, 5, 6], [1, 2, 3, 4, 5, 6])
        if (allocated(error)) return

    contains
        subroutine check_to_array(s, e, expected)
            integer, intent(in) :: s(3), e(3), expected(6)
            type(t_Block) :: blk

            integer :: tarr(6) = 0

            blk = new_Block(s, e)
            call blk%to_array(tarr)
            call check(error, tarr, expected)
            if (allocated(error)) return
        end subroutine
    end subroutine

    subroutine test_from_array(error)
        type(error_type), allocatable, intent(out) :: error

        call check_from_array([0, 0, 0, 0, 0, 0], [0, 0, 0], [0, 0, 0])
        if (allocated(error)) return

        call check_from_array([1, 2, 3, 4, 5, 6], [1, 2, 3], [4, 5, 6])
        if (allocated(error)) return

        call check_from_array([6, 5, 4, 3, 2, 1], [6, 5, 4], [3, 2, 1])
        if (allocated(error)) return

    contains
        subroutine check_from_array(arr, s, e)
            integer, intent(in) :: arr(6), s(3), e(3)
            type(t_Block) :: blk

            integer :: tarr(6) = 0

            call blk%from_array(arr)

            call check(error, blk%start, s)
            if (allocated(error)) return

            call check(error, blk%end, e)
            if (allocated(error)) return
        end subroutine
    end subroutine

end module
