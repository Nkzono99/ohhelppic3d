module m_mpi_fft_solver
    use m_block
    implicit none

    private
    public t_MPIFFTSolver3d

    type, abstract :: t_MPIFFTSolver3d
        integer :: boundary_types(3)
        type(t_Block) :: local_block
        type(t_Block) :: logical_block
        type(t_Block) :: global_block
    contains
        procedure(mpiFFTSolver3d_forward), deferred :: forward
        procedure(mpiFFTSolver3d_backward), deferred :: backward
    end type

    interface
        subroutine mpiFFTSolver3d_forward(self, in, out)
            import t_MPIFFTSolver3d
            class(t_MPIFFTSolver3d), intent(inout) :: self
            double precision, intent(in) :: in(self%local_block%start(1):self%local_block%end(1), &
                                               self%local_block%start(2):self%local_block%end(2), &
                                               self%local_block%start(3):self%local_block%end(3))
            double precision, intent(out) :: out(self%local_block%start(1):self%local_block%end(1), &
                                                 self%local_block%start(2):self%local_block%end(2), &
                                                 self%local_block%start(3):self%local_block%end(3))
        end subroutine

        subroutine mpiFFTSolver3d_backward(self, in, out)
            import t_MPIFFTSolver3d
            class(t_MPIFFTSolver3d), intent(inout) :: self
            double precision, intent(in) :: in(self%local_block%start(1):self%local_block%end(1), &
                                               self%local_block%start(2):self%local_block%end(2), &
                                               self%local_block%start(3):self%local_block%end(3))
            double precision, intent(out) :: out(self%local_block%start(1):self%local_block%end(1), &
                                                 self%local_block%start(2):self%local_block%end(2), &
                                                 self%local_block%start(3):self%local_block%end(3))
        end subroutine
    end interface

end module
