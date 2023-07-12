module m_poisson_field_solver
    use m_mpi_fft_solver
    use m_science_constants, only: pi
    use m_get_default, only: get_default
    use m_block
    use m_field_boundary_type
    use m_field_solver
    use m_ohfield
    use m_ohhelp
    use m_string_holder
    use m_mpi_fft_solver_factory
    use m_poisson_solver3d
    implicit none

    !> 3d poisson equation solver.
    !>
    !> Poisson equation:
    !>     ∂^2p/∂^2 + ∂^2p/∂y^2 + ∂^2p/∂z^2 = f(x, y, z)
    type, extends(t_FieldSolver) :: t_PoissonFieldSolver3d
        class(t_MPIFFTSolver3d), pointer, private :: fft3d
        type(t_PoissonSolver3d), private :: poisson_solver3d
        double precision, allocatable, private :: modified_wave_number(:, :, :)
        double precision, private :: boundary_condition_terms(2, 3)
        type(t_Block) :: local_block
        type(t_Block) :: global_block
    contains
        procedure :: solve => poissonSolver3d_solve
    end type

    private
    public new_PoissonFieldSolver3d

contains

    function new_PoissonFieldSolver3d(local_block, global_block, &
                                      fft_solver_name, &
                                      boundary_types, &
                                      boundary_values, &
                                      myid, &
                                      nprocs, &
                                      comm, &
                                      tag) result(obj)
        type(t_Block), intent(in) :: local_block
        type(t_Block), intent(in) :: global_block
        character(*), intent(in) :: fft_solver_name
        type(t_StringHolder), intent(in) :: boundary_types(3)
        double precision, intent(in), optional :: boundary_values(2, 3)
        integer, intent(in) :: myid
        integer, intent(in) :: nprocs
        integer, intent(in) :: comm
        integer, intent(in) :: tag
        type(t_PoissonFieldSolver3d) :: obj

        integer :: start(3), end(3)
        integer :: kx, ky, kz

        block
            integer :: fft_boundary_types(3)

            integer :: i

            do i = 1, 3
                select case (boundary_types(i)%string)
                case ('periodic')
                    fft_boundary_types(i) = Field_BoundaryType_Periodic
                case ('dirichlet')
                    fft_boundary_types(i) = Field_BoundaryType_Dirichlet
                case ('neumman')
                    fft_boundary_types(i) = Field_BoundaryType_Neumann
                case ('dirichlet-neumman')
                    fft_boundary_types(i) = Field_BoundaryType_Dirichlet_Neumann
                case ('neumman-dirichlet')
                    fft_boundary_types(i) = Field_BoundaryType_Neumann_Dirichlet
                end select
            end do

            allocate (obj%fft3d, source=create_mpi_fft_solver(fft_solver_name, &
                                                              fft_boundary_types, &
                                                              local_block, &
                                                              global_block, &
                                                              myid, nprocs, &
                                                              comm, tag=10))
        end block

        obj%local_block = local_block
        obj%global_block = global_block

        obj%poisson_solver3d = new_PoissonSolver3d(local_block, &
                                                   global_block, &
                                                   obj%fft3d)
    end function

    subroutine poissonSolver3d_solve(self, rho, aj, eb, phi, ohhelp)
        class(t_PoissonFieldSolver3d), intent(in) :: self
        class(t_OhField), intent(in) :: rho
        class(t_OhField), intent(in) :: aj
        class(t_OhField), intent(inout) :: eb
        class(t_OhField), intent(inout) :: phi
        class(t_OhHelp), intent(inout) :: ohhelp

        eb%values = 0

        block
            integer :: local_start(3), local_end(3)

            local_start = phi%to_local_index(self%local_block%start, 1)
            local_end = phi%to_local_index(self%local_block%end, 1)

            call self%poisson_solver3d%solve(rho%values(1, &
                                                        local_start(1):local_end(1), &
                                                        local_start(2):local_end(2), &
                                                        local_start(3):local_end(3), &
                                                        1), &
                                             phi%values(1, &
                                                        local_start(1):local_end(1), &
                                                        local_start(2):local_end(2), &
                                                        local_start(3):local_end(3), &
                                                        1))

            call ohhelp%broadcast_field(phi)
            call ohhelp%exchange_borders(phi)
        end block

        ! NOTE: ebをbroadcastせず、primary/secondary両方で解いているのは通信回数を減らすことを意図している
        ! OPTIMIZE: しかし、実行時間の計測・比較はしていないため、どちらを採用するかは要検討
        block
            integer :: ps
            integer :: i, j, k
            integer :: xl, yl, zl
            integer :: xu, yu, zu

            do ps = 1, 2
                xl = phi%subdomain_range(1, 1, ps) - 1
                yl = phi%subdomain_range(1, 2, ps) - 1
                zl = phi%subdomain_range(1, 3, ps) - 1
                xu = phi%subdomain_range(2, 1, ps) - phi%subdomain_range(1, 1, ps)
                yu = phi%subdomain_range(2, 2, ps) - phi%subdomain_range(1, 2, ps)
                zu = phi%subdomain_range(2, 3, ps) - phi%subdomain_range(1, 3, ps)

                eb%values(1, xl:xu, yl:yu, zl:zu, ps) = &
                    phi%values(1, xl:xu, yl:yu, zl:zu, ps) - phi%values(1, xl + 1:xu + 1, yl:yu, zl:zu, ps)
                eb%values(2, xl:xu, yl:yu, zl:zu, ps) = &
                    phi%values(2, xl:xu, yl:yu, zl:zu, ps) - phi%values(2, xl:xu, yl + 1:yu + 1, zl:zu, ps)
                eb%values(3, xl:xu, yl:yu, zl:zu, ps) = &
                    phi%values(3, xl:xu, yl:yu, zl:zu, ps) - phi%values(3, xl:xu, yl:yu, zl + 1:zu + 1, ps)
            end do

            call ohhelp%exchange_borders(eb)
        end block

        ! Avoiding self-force.
        block
            integer :: ps
            integer :: i, j, k
            integer :: xl, yl, zl
            integer :: xu, yu, zu

            do ps = 1, 2
                xl = phi%subdomain_range(1, 1, ps) - 1
                yl = phi%subdomain_range(1, 2, ps) - 1
                zl = phi%subdomain_range(1, 3, ps) - 1
                xu = phi%subdomain_range(2, 1, ps) - phi%subdomain_range(1, 1, ps)
                yu = phi%subdomain_range(2, 2, ps) - phi%subdomain_range(1, 2, ps)
                zu = phi%subdomain_range(2, 3, ps) - phi%subdomain_range(1, 3, ps)

                eb%values(1, xl:xu, yl:yu, zl:zu, ps) = &
                    0.5d0*(eb%values(1, xl - 1:xu - 1, yl:yu, zl:zu, ps) + eb%values(1, xl:xu, yl:yu, zl:zu, ps))
                eb%values(2, xl:xu, yl:yu, zl:zu, ps) = &
                    0.5d0*(eb%values(2, xl:xu, yl - 1:yu - 1, zl:zu, ps) + eb%values(2, xl:xu, yl:yu, zl:zu, ps))
                eb%values(3, xl:xu, yl:yu, zl:zu, ps) = &
                    0.5d0*(eb%values(3, xl:xu, yl:yu, zl - 1:zu - 1, ps) + eb%values(3, xl:xu, yl:yu, zl:zu, ps))
            end do

            call ohhelp%exchange_borders(eb)
        end block
    end subroutine

end module
