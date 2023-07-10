module test_pcg_generator
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use m_pcg_generator
    use m_random_generator
    use m_check_array, only: check
    implicit none

    private
    public collect_pcg_generator

contains

    subroutine collect_pcg_generator(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("pcg_generator.rand", test_rand), &
                    new_unittest('pcg_generator.advance', test_advance) &
                    ]
    end subroutine

    subroutine test_rand(error)
        type(error_type), allocatable, intent(out) :: error

        class(t_RandomGenerator), allocatable :: random_generator

        random_generator = new_PcgGenerator([int(42, kind=8), int(50, kind=8)])

        print *, 'rand(): [0, 1.0)'
        block
            integer :: i
            do i = 1, 10
                print *, random_generator%rand()
            end do
        end block

        block
            integer :: i
            double precision :: sum = 0d0

            do i = 1, 1000
                sum = sum + random_generator%rand()
            end do

            print *, 'mean(rand() for i in range(1000)) = ', sum/1000
        end block
    end subroutine

    subroutine test_advance(error)
        type(error_type), allocatable, intent(out) :: error

        class(t_RandomGenerator), allocatable :: random_generator

        random_generator = new_PcgGenerator([int(42, kind=8), int(50, kind=8)])

        print *, 'rand(): [0, 1.0)'
        block
            integer :: i
            do i = 1, 5
                print *, random_generator%rand()
            end do
        end block

        print *, 'call advance(-5)'
        call random_generator%advance(-5)

        print *, 'rand(): [0, 1.0)'
        block
            integer :: i
            do i = 1, 5
                print *, random_generator%rand()
            end do
        end block

        print *, 'call advance(-5)'
        call random_generator%advance(-5)
        print *, 'call advance(1000000000)'
        call random_generator%advance(1000000000)

        print *, 'rand(): [0, 1.0)'
        block
            integer :: i
            do i = 1, 5
                print *, random_generator%rand()
            end do
        end block
    end subroutine
end module
