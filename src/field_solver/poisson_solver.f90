module m_poisson_solver
    use m_mpi_fft_solver
    use m_science_constants, only: pi
    use m_get_default, only: get_default
    use m_block
    use m_field_boundary_type
    use m_field_solver
    use m_ohfield
    use m_ohhelp
    implicit none

    !> 3d poisson equation solver.
    !>
    !> Poisson equation:
    !>     ∂^2p/∂^2 + ∂^2p/∂y^2 + ∂^2p/∂z^2 = f(x, y, z)
    type, extends(t_FieldSolver) :: t_PoissonSolver3d
        class(t_MPIFFTSolver3d), pointer, private :: fft3d
        double precision, allocatable, private :: modified_wave_number(:, :, :)
        double precision, private :: boundary_condition_terms(2, 3)
        type(t_Block) :: local_block
        type(t_Block) :: global_block
    contains
        procedure :: solve => poissonSolver3d_solve
    end type

    private
    public new_PoissonSolver3d

contains

    function new_PoissonSolver3d(local_block, global_block, fft3d, boundary_values) result(obj)
        type(t_Block), intent(in) :: local_block
        type(t_Block), intent(in) :: global_block
        class(t_MPIFFTSolver3d), pointer, intent(in) :: fft3d
        double precision, intent(in), optional :: boundary_values(2, 3)
        type(t_PoissonSolver3d) :: obj

        integer :: start(3), end(3)
        integer :: kx, ky, kz

        obj%fft3d => fft3d

        obj%local_block = local_block
        obj%global_block = global_block

        if (present(boundary_values)) then
            obj%boundary_condition_terms(:, 1) = &
                calc_boundary_term(boundary_values(:, 1), obj%fft3d%boundary_types(1))
            obj%boundary_condition_terms(:, 2) = &
                calc_boundary_term(boundary_values(:, 2), obj%fft3d%boundary_types(2))
            obj%boundary_condition_terms(:, 3) = &
                calc_boundary_term(boundary_values(:, 3), obj%fft3d%boundary_types(3))
        else
            obj%boundary_condition_terms(:, :) = 0.0d0
        end if

        allocate (obj%modified_wave_number(local_block%start(1):local_block%end(1), &
                                           local_block%sizes(2):local_block%end(2), &
                                           local_block%sizes(3):local_block%end(3)))

        start(:) = local_block%start(:) - global_block%start(:)
        end(:) = local_block%end(:) - global_block%end(:)
        do concurrent(kx=start(1):end(1), ky=start(2):end(2), kz=start(3):end(3))
            block
                double precision :: wx, wy, wz
                double precision :: wn

                wx = calc_wave_number(kx, global_block%sizes(1), fft3d%boundary_types(1))
                wy = calc_wave_number(ky, global_block%sizes(2), fft3d%boundary_types(2))
                wz = calc_wave_number(kz, global_block%sizes(3), fft3d%boundary_types(3))
                wn = wx + wy + wz

                obj%modified_wave_number(kx + global_block%start(1), &
                                         ky + global_block%start(2), &
                                         kz + global_block%start(3)) = wn
            end block
        end do
    end function

    pure function calc_boundary_term(boundary_values, boundary_type) result(terms)
        double precision, intent(in) :: boundary_values(2)
        integer, intent(in) :: boundary_type
        double precision :: terms(2)

        select case (boundary_type)
        case (Field_BoundaryType_Periodic)
            terms(:) = boundary_values(:)

        case (Field_BoundaryType_Dirichlet)
            terms(:) = [-boundary_values(1), -boundary_values(2)]

        case (Field_BoundaryType_Neumann)
            terms(:) = [2.0d0*boundary_values(1), -2.0d0*boundary_values(2)]

        case (Field_BoundaryType_Dirichlet_Neumann)
            terms(:) = [-boundary_values(1), -2.0d0*boundary_values(2)]

        case (Field_BoundaryType_Neumann_Dirichlet)
            terms(:) = [-boundary_values(1), -2.0d0*boundary_values(2)]
        end select
    end function

    pure function calc_wave_number(k, n, boundary_type) result(wn)
        integer, intent(in) :: k
        integer, intent(in) :: n
        integer, intent(in) :: boundary_type
        double precision :: wn

        select case (boundary_type)
        case (Field_BoundaryType_Periodic)
            if (k <= int(n/2)) then
                wn = 2.0d0*sin(PI*k/dble(n))
            else
                wn = 2.0d0*sin(PI*(k - int(n/2))/dble(n))
            end if

        case (Field_BoundaryType_Dirichlet)
            wn = 2.0d0*(cos(PI*(k + 1)/dble(n + 1)) - 1.0d0)

        case (Field_BoundaryType_Neumann)
            wn = 2.0d0*(cos(PI*k/dble(n)) - 1.0d0)

        case (Field_BoundaryType_Dirichlet_Neumann)
            wn = 2.0d0*(cos(PI*(k + 0.5d0)/dble(n + 1)) - 1.0d0)

        case (Field_BoundaryType_Neumann_Dirichlet)
            wn = 2.0d0*(cos(PI*(k + 0.5d0)/dble(n + 1)) - 1.0d0)

        end select
    end function

    subroutine poissonSolver3d_solve(self, rho, aj, eb, phi, ohhelp)
        class(t_PoissonSolver3d), intent(in) :: self
        class(t_OhField), intent(in) :: rho
        class(t_OhField), intent(in) :: aj
        class(t_OhField), intent(inout) :: eb
        class(t_OhField), intent(inout) :: phi
        class(t_OhHelp), intent(inout) :: ohhelp

        eb%values = 0

        phi%values(1, :, :, :, 1) = rho%values(1, :, :, :, 1)

        block
            integer :: axis
            integer :: local_start(3), local_end(3)

            local_start = phi%to_local_index(self%global_block%start, 1)
            local_end = phi%to_local_index(self%global_block%end, 1)

            do axis = 1, 3
                if (self%local_block%start(axis) == self%global_block%start(axis)) then
                    phi%values(1, local_start(axis), :, :, 1) = &
                        phi%values(1, local_start(axis), :, :, 1) &
                        + self%boundary_condition_terms(1, axis)
                end if

                if (self%local_block%end(axis) == self%global_block%end(axis)) then
                    phi%values(1, local_end(axis), :, :, 1) = &
                        phi%values(1, local_end(axis), :, :, 1) &
                        + self%boundary_condition_terms(2, axis)
                end if
            end do
        end block

        block
            integer :: local_start(3), local_end(3)

            local_start = phi%to_local_index(self%local_block%start, 1)
            local_end = phi%to_local_index(self%local_block%end, 1)

            call self%fft3d%forward(phi%values(1, &
                                               local_start(1):local_end(1), &
                                               local_start(2):local_end(2), &
                                               local_start(3):local_end(3), &
                                               1), &
                                    phi%values(1, &
                                               local_start(1):local_end(1), &
                                               local_start(2):local_end(2), &
                                               local_start(3):local_end(3), &
                                               1))

            phi%values(1, &
                       local_start(1):local_end(1), &
                       local_start(2):local_end(2), &
                       local_start(3):local_end(3), &
                       1) = &
                phi%values(1, &
                           local_start(1):local_end(1), &
                           local_start(2):local_end(2), &
                           local_start(3):local_end(3), &
                           1) &
                /self%modified_wave_number(self%local_block%start(1):self%local_block%end(1), &
                                           self%local_block%start(2):self%local_block%end(2), &
                                           self%local_block%start(3):self%local_block%end(3))

            if (all(self%local_block%start == self%global_block%start)) then
                phi%values(1, &
                           local_start(1), &
                           local_start(2), &
                           local_start(3), &
                           1) = 0d0
            end if

            call self%fft3d%backward(phi%values(1, &
                                                local_start(1):local_end(1), &
                                                local_start(2):local_end(2), &
                                                local_start(3):local_end(3), &
                                                1), &
                                     phi%values(1, &
                                                local_start(1):local_end(1), &
                                                local_start(2):local_end(2), &
                                                local_start(3):local_end(3), &
                                                1))
        end block

        call ohhelp%broadcast_field(phi)

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
        end block
    end subroutine

end module
