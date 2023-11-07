module m_mpi_fft_solver_factory
    use m_block
    use m_mpi_fft_solver
    use m_mpi_fftw3_solver
    implicit none

    private
    public create_mpi_fft_solver

contains

    function create_mpi_fft_solver(fft_solver_name, &
                                   fft_boundary_types, &
                                   local_block, &
                                   global_block, &
                                   myid, &
                                   nprocs, &
                                   comm, &
                                   tag) result(obj)
        type(t_Block), intent(in) :: local_block
        type(t_Block), intent(in) :: global_block
        integer, intent(in) :: fft_boundary_types(3)
        character(*), intent(in) :: fft_solver_name
        integer, intent(in) :: myid
        integer, intent(in) :: nprocs
        integer, intent(in) :: comm
        integer, intent(in) :: tag
        class(t_MPIFFTSolver3d), allocatable :: obj

        select case (fft_solver_name)
        case ('fftw3')
            obj = new_MPIFFTW3Solver3d(fft_boundary_types, &
                                       local_block, &
                                       global_block, &
                                       myid, &
                                       nprocs, &
                                       comm, &
                                       tag=tag)
        end select
    end function

end module
