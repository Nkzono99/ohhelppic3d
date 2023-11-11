module test_poisson_solver
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use m_poisson_solver3d
    use m_check_array, only: check
    use m_mpi_fft_solver
    use m_mpi_fftw3_solver
    use m_field_boundary_type
    use m_block
    use mpi
    use m_domain
    use m_hdf5
    implicit none

    private
    public collect_poisson_solver

contains

    subroutine collect_poisson_solver(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("new_PoissonSolver", test_new_PoissonSolver) &
                    ]
    end subroutine

    subroutine test_new_PoissonSolver(error)
        type(error_type), allocatable, intent(out) :: error
        type(t_PoissonSolver3d) :: obj
        class(t_MPIFFTSolver3d), allocatable, target :: fft3d

        integer :: n(3), nnodes(3)
        integer :: nprocs, myid
        integer :: ierr

        call mpi_init(ierr)

        call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
        call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)

        n = [32, 32, 32]
        nnodes = [2, 1, 2]

        block
            integer :: boundary_types(3)
            double precision, allocatable :: rho(:, :, :), phi(:, :, :)

            type(t_Block) :: local_block, global_block
            integer :: start(3), end(3), i

            integer :: myid_xyz(3)

            myid_xyz = [ &
                       mod(myid, nnodes(1)), &
                       mod(myid/nnodes(1), nnodes(2)), &
                       mod(myid/nnodes(1)/nnodes(2), nnodes(3)) &
                       ]

            start(:) = myid_xyz*(n/nnodes)
            end(:) = (myid_xyz + 1)*(n/nnodes)

            do i = 1, 3
                if (myid_xyz(i) == nnodes(i) - 1) end(i) = n(i)
            end do

            local_block = new_Block(start, end)
            global_block = new_Block([0, 0, 0], n)

            boundary_types = [ &
                             Field_BoundaryType_Periodic, &
                             Field_BoundaryType_Periodic, &
                             Field_BoundaryType_Periodic &
                             ]
            fft3d = new_MPIFFTW3Solver3d(boundary_types, &
                                         local_block, global_block, &
                                         myid, nprocs, MPI_COMM_WORLD, 10)

            obj = new_PoissonSolver3d(local_block, global_block, fft3d)

            allocate (rho( &
                      0:local_block%sizes(1)-1, &
                      0:local_block%sizes(2)-1, &
                      0:local_block%sizes(3)-1 &
                      ))
            allocate (phi( &
                      0:local_block%sizes(1)-1, &
                      0:local_block%sizes(2)-1, &
                      0:local_block%sizes(3)-1 &
                      ))

            rho(:, :, :) = 0
            block
                integer :: p(3)

                p = [8, 10, 10]
                if (local_block%start(1) <= p(1) .and. p(1) <= local_block%end(1) &
                    .and. local_block%start(2) <= p(2) .and. p(2) <= local_block%end(2) &
                    .and. local_block%start(3) <= p(3) .and. p(3) <= local_block%end(3)) then
                    rho( &
                        p(1) - local_block%start(1), &
                        p(2) - local_block%start(2), &
                        p(3) - local_block%start(3) &
                        ) = -1.0d0
                end if
            end block

            call obj%solve(rho, phi)

            block
                type(t_HDF5File) :: h5
                type(t_HDF5Group) :: group
                type(t_SubDomain3d) :: subdomain_info
                integer :: status

                subdomain_info = new_SubDomain3d( &
                                 int(local_block%sizes, kind=8), &
                                 int(local_block%start, kind=8), &
                                 int(global_block%sizes, kind=8), &
                                 int(global_block%start, kind=8) &
                                 )

                call hdf5_initialize(status)

                h5 = new_HDF5File('phisp00_0000.h5', 'w', MPI_COMM_WORLD)
                group = h5%create_group('phisp')

                call group%write_dataset('0', &
                                         phi, subdomain_info, status)

                call group%close(status)
                call h5%close(status)

                call hdf5_finalize(status)
            end block
        end block

    end subroutine
end module
