module m_mpi_fft_solver
    implicit none

    private
    public t_MPIFFTSolver3d

    type :: t_MPIFFTSolver3d
        integer :: boundary_types(3)
    contains
        procedure :: forward => mpiFFTSolver3d_forward
        procedure :: backward => mpiFFTSolver3d_backward
    end type

    interface
        subroutine mpiFFTSolver3d_forward(self, in, out)
            import t_MPIFFTSolver3d
            class(t_MPIFFTSolver3d), intent(inout) :: self
            double precision :: in(:, :, :)
            double precision :: out(:, :, :)
        end subroutine

        subroutine mpiFFTSolver3d_backward(self, in, out)
            import t_MPIFFTSolver3d
            class(t_MPIFFTSolver3d), intent(inout) :: self
            double precision :: in(:, :, :)
            double precision :: out(:, :, :)
        end subroutine
    end interface

end module
