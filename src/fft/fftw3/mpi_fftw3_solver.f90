module m_mpi_fftw_solver
    use m_mpi_fft_solver
    use m_fftw3_mpi
    use m_fft_boundary_type
    use m_block
    use m_mpi_block_rebaser
    use mpi
    implicit none

    private
    public t_MPIFFTWSolver3d
    public new_MPIFFTWSolver3d

    type, extends(t_MPIFFTSolver3d) :: t_MPIFFTWSolver3d
        type(C_PTR) :: forward_plan
        type(C_PTR) :: backward_plan
        double precision :: normalizers(3)

        type(t_Block) :: require_block

        integer :: myid
        integer :: nprocs

        integer :: communicator
        integer :: tag

        type(C_PTR) :: ptr_fft_array
        real(C_DOUBLE), pointer :: fft_array(:, :, :)

        type(t_MPIBlockRebaser) :: in2fftarray
        type(t_MPIBlockRebaser) :: fftarray2out
    contains
        procedure :: forward => mpiFFTWSolver_forward
        procedure :: backward => mpiFFTWSolver_backward
        procedure :: cleanup => mpiFFTWSolver_cleanup
    end type

contains

    subroutine convert_fft_type(n, boundary_type, forward_fft_type, backward_fft_type, normalizer)
        integer(C_INT), intent(in) :: n
        integer, intent(in) :: boundary_type
        integer(C_FFTW_R2R_KIND), intent(out) :: forward_fft_type
        integer(C_FFTW_R2R_KIND), intent(out) :: backward_fft_type
        double precision, intent(out) :: normalizer

        select case (boundary_type)
        case (FFT_BoundaryType_Periodic)
            forward_fft_type = FFTW_R2HC
            backward_fft_type = FFTW_HC2R
            normalizer = n

        case (FFT_BoundaryType_Dirichlet)
            forward_fft_type = FFTW_RODFT00
            backward_fft_type = FFTW_RODFT00
            normalizer = 2.0d0*(n + 1.0d0)

        case (FFT_BoundaryType_Neumann)
            forward_fft_type = FFTW_REDFT00
            backward_fft_type = FFTW_REDFT00
            normalizer = 2.0d0*(n - 1.0d0)

        case (FFT_BoundaryType_Dirichlet_Neumann)
            forward_fft_type = FFTW_RODFT01
            backward_fft_type = FFTW_RODFT10
            normalizer = 2.0d0*n

        case (FFT_BoundaryType_Neumann_Dirichlet)
            forward_fft_type = FFTW_REDFT01
            backward_fft_type = FFTW_REDFT10
            normalizer = 2.0d0*n
        end select
    end subroutine

    function new_MPIFFTWSolver3d(fft_boundary_types, &
                                 local_block, &
                                 global_block, &
                                 myid, nprocs, &
                                 communicator, tag) result(obj)
        integer, intent(in) :: fft_boundary_types(3)
        !> Block that the local process has.
        !> global index 座標系は1:nx
        type(t_Block), intent(in) :: local_block
        type(t_Block), intent(in) :: global_block
        integer, intent(in) :: myid
        integer, intent(in) :: nprocs
        integer, intent(in) :: communicator
        integer, intent(in) :: tag
        type(t_MPIFFTWSolver3d) :: obj

        integer :: nx, ny, nz
        integer(kind=C_INTPTR_T) :: iptr_nx, iptr_ny, iptr_nz
        integer(C_INTPTR_T) :: local_nz, local_z_start

        obj%boundary_types(:) = fft_boundary_types(:)
        obj%local_block = local_block
        obj%global_block = global_block
        obj%myid = myid
        obj%nprocs = nprocs
        obj%communicator = communicator
        obj%tag = tag

        nx = global_block%sizes(1)
        ny = global_block%sizes(2)
        nz = global_block%sizes(3)

        iptr_nx = nx
        iptr_ny = ny
        iptr_nz = nz

        call fftw_mpi_init()

        block ! Allocate array used in fft.
            integer(C_INTPTR_T) :: alloc_local

            alloc_local = fftw_mpi_local_size_3d(iptr_nz, &
                                                 iptr_ny, &
                                                 iptr_nx, &
                                                 communicator, &
                                                 local_nz, local_z_start)
            obj%ptr_fft_array = fftw_alloc_real(alloc_local)
            call c_f_pointer(obj%ptr_fft_array, obj%fft_array, &
                             [iptr_nx, iptr_ny, local_nz])
        end block

        block ! Create fftw plan.
            integer(C_FFTW_R2R_KIND) :: forward_fft_type_x
            integer(C_FFTW_R2R_KIND) :: forward_fft_type_y
            integer(C_FFTW_R2R_KIND) :: forward_fft_type_z
            integer(C_FFTW_R2R_KIND) :: backward_fft_type_x
            integer(C_FFTW_R2R_KIND) :: backward_fft_type_y
            integer(C_FFTW_R2R_KIND) :: backward_fft_type_z
            double precision :: normalizer_x
            double precision :: normalizer_y
            double precision :: normalizer_z

            call convert_fft_type(nx, fft_boundary_types(1), &
                                  forward_fft_type_x, backward_fft_type_x, normalizer_x)
            call convert_fft_type(ny, fft_boundary_types(1), &
                                  forward_fft_type_y, backward_fft_type_y, normalizer_y)
            call convert_fft_type(nz, fft_boundary_types(1), &
                                  forward_fft_type_z, backward_fft_type_z, normalizer_z)

            obj%normalizers = [normalizer_x, normalizer_y, normalizer_z]

            obj%forward_plan = fftw_mpi_plan_r2r_3d(iptr_nz, iptr_ny, iptr_nx, &
                                                    obj%fft_array, obj%fft_array, &
                                                    obj%communicator, &
                                                    forward_fft_type_z, &
                                                    forward_fft_type_y, &
                                                    forward_fft_type_x, &
                                                    FFTW_MEASURE)

            obj%backward_plan = fftw_mpi_plan_r2r_3d(iptr_nz, iptr_ny, iptr_nx, &
                                                     obj%fft_array, obj%fft_array, &
                                                     obj%communicator, &
                                                     backward_fft_type_z, &
                                                     backward_fft_type_y, &
                                                     backward_fft_type_x, &
                                                     FFTW_MEASURE)
        end block

        block ! Set the coordinate range (block) to be used in fftw.
            integer :: zs, ze

            ! local_z from fftw start at 0 to nz-1
            zs = local_z_start + global_block%start(3)
            ze = local_z_start + global_block%start(3) + local_nz - 1
            obj%require_block = new_Block([global_block%start(1), global_block%start(2), zs], &
                                          [global_block%end(1), global_block%end(2), ze])
        end block

        block ! Create a converter between the coordinate range used in fft and the coordinate range for subdomain.
            integer :: pids(nprocs)
            integer :: i

            do i = 1, nprocs
                pids(i) = i - 1
            end do

            obj%in2fftarray = new_MPIBlockRebaser(obj%local_block, obj%require_block, pids, myid + 1, obj%communicator)
            obj%fftarray2out = new_MPIBlockRebaser(obj%require_block, obj%local_block, pids, myid + 1, obj%communicator)
        end block
    end function

    subroutine mpiFFTWSolver_forward(self, in, out)
        class(t_MPIFFTWSolver3d), intent(inout) :: self
        double precision, intent(inout) :: in(self%local_block%start(1):self%local_block%end(1), &
                                              self%local_block%start(2):self%local_block%end(2), &
                                              self%local_block%start(3):self%local_block%end(3))
        double precision, intent(inout) :: out(self%local_block%start(1):self%local_block%end(1), &
                                               self%local_block%start(2):self%local_block%end(2), &
                                               self%local_block%start(3):self%local_block%end(3))

        call self%in2fftarray%rebase(in, self%fft_array(:, :, :), tag=self%tag)
        call fftw_mpi_execute_r2r(self%forward_plan, self%fft_array, self%fft_array)
        call self%fftarray2out%rebase(self%fft_array(:, :, :), out, tag=self%tag)
    end subroutine

    subroutine mpiFFTWSolver_backward(self, in, out)
        class(t_MPIFFTWSolver3d), intent(inout) :: self
        double precision, intent(inout) :: in(self%local_block%start(1):self%local_block%end(1), &
                                              self%local_block%start(2):self%local_block%end(2), &
                                              self%local_block%start(3):self%local_block%end(3))
        double precision, intent(inout) :: out(self%local_block%start(1):self%local_block%end(1), &
                                               self%local_block%start(2):self%local_block%end(2), &
                                               self%local_block%start(3):self%local_block%end(3))

        call self%in2fftarray%rebase(in, self%fft_array, tag=self%tag)

        call fftw_mpi_execute_r2r(self%backward_plan, self%fft_array, self%fft_array)

        call self%fftarray2out%rebase(self%fft_array, out, tag=self%tag)

        out = out/product(self%normalizers)
    end subroutine

    subroutine mpiFFTWSolver_cleanup(self)
        class(t_MPIFFTWSolver3d), intent(inout) :: self

        call fftw_free(self%ptr_fft_array)
        call fftw_mpi_cleanup()
    end subroutine

end module
