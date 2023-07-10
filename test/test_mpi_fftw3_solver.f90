#define WITH_XDP 1
module test_mpi_fftw3_solver
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use m_mpi_fftw3_solver
    use m_mpi_fft_solver
    use m_field_boundary_type
    use mpi
    use m_block
    use m_mpi_block
    use m_check_array, only: check
    use m_science_constants
    use m_str
    implicit none

    private
    public collect_mpi_fftw3_solver

    integer :: nprocs, myid

contains

    subroutine collect_mpi_fftw3_solver(testsuite)
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
                    new_unittest("mpiFftw3Solver.forward [myid="//str(myid)//"]", test_forward) &
                    ]
    end subroutine

    subroutine test_forward(error)
        type(error_type), allocatable, intent(out) :: error
        class(t_MPIFFTSolver3d), allocatable :: fft_solver
        type(t_Block) :: local_block, global_block

        integer :: boundary_types(3)
        integer :: n(3), nnodes(3), ids(3)

        double precision, allocatable :: in(:, :, :), out(:, :, :)
        double precision, allocatable :: expected(:)

        ! Test impulse response.
        n = [1, 1, 63]
        nnodes = [1, 1, 4]
        boundary_types(:) = Field_BoundaryType_Periodic
        global_block = new_Block([0, 0, 0], n)

        call create_local_block(global_block, nnodes, myid, local_block, ids)

        fft_solver = new_MPIFFTW3Solver3d(boundary_types, &
                                          local_block, global_block, &
                                          myid, nprocs, MPI_COMM_WORLD, 100)

        allocate (in(local_block%sizes(1), local_block%sizes(2), local_block%sizes(3)))
        allocate (out(local_block%sizes(1), local_block%sizes(2), local_block%sizes(3)))

        in = 0
        if (myid == 0) then
            in(1, 1, 2) = 1
        end if

        call fft_solver%forward(in, out)

        allocate (expected(local_block%sizes(3)))

        block
            double precision :: x
            integer :: i

            do i = local_block%start(3), local_block%end(3)
                x = i

                ! FFTW_R2HC: r0, r1, r2, ..., r_{n/2}, i_{(n+1)/2-1}, ..., i2, i1
                if (i <= (global_block%sizes(3) + 1)/2 - 1) then
                    expected(i - local_block%start(3) + 1) = cos(x/(global_block%sizes(3) - 1)*2*pi)
                else
                    expected(i - local_block%start(3) + 1) = sin(x/(global_block%sizes(3) - 1)*2*pi)
                end if
            end do
        end block

        block
            double precision :: diff

            diff = sum(abs(out(1, 1, :) - expected(:)))
            call check(error, diff, 0d0, thr=1d-8)
        end block
    end subroutine
end module
