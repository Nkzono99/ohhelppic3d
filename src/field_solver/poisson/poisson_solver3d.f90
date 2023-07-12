module m_poisson_solver3d
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
    implicit none

    !> 3d poisson equation solver.
    !>
    !> Poisson equation:
    !>     ∂^2p/∂^2 + ∂^2p/∂y^2 + ∂^2p/∂z^2 = f(x, y, z)
    type :: t_PoissonSolver3d
        class(t_MPIFFTSolver3d), pointer, private :: fft3d
        double precision, allocatable, private :: modified_wave_number(:, :, :)
        double precision, private :: boundary_condition_terms(2, 3)
        type(t_Block) :: local_block
        type(t_Block) :: global_block
    contains
        procedure :: solve => poissonSolver3d_solve
    end type

    private
    public t_PoissonSolver3d
    public new_PoissonSolver3d

contains

    function new_PoissonSolver3d(local_block, global_block, &
                                 fft3d, &
                                 boundary_values) result(obj)
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
                                           local_block%start(2):local_block%end(2), &
                                           local_block%start(3):local_block%end(3)))

        start(:) = local_block%start(:) - global_block%start(:)
        end(:) = local_block%end(:) - global_block%start(:)
        do concurrent(kx=start(1):end(1), ky=start(2):end(2), kz=start(3):end(3))
            block
                double precision :: wx, wy, wz
                double precision :: wn

                wx = calc_wave_number(kx, global_block%sizes(1) - 1, obj%fft3d%boundary_types(1))
                wy = calc_wave_number(ky, global_block%sizes(2) - 1, obj%fft3d%boundary_types(2))
                wz = calc_wave_number(kz, global_block%sizes(3) - 1, obj%fft3d%boundary_types(3))

                ! TODO: 二乗にすると正しそうな結果となったが、理論をもう一度確認すること
                wn = wx*wx + wy*wy + wz*wz

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
            ! REVIEW: wave numberの計算が本当に正しいか疑わしいため、要再確認
            if (k <= int(n/2)) then
                wn = 2.0d0*sin(PI*k/dble(n))
            else
                wn = 2.0d0*sin(PI*(n - k)/dble(n))
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

    subroutine poissonSolver3d_solve(self, rho, phi)
        class(t_PoissonSolver3d), intent(in) :: self
        double precision, intent(in) :: rho( &
            self%local_block%start(1):self%local_block%end(1), &
            self%local_block%start(2):self%local_block%end(2), &
            self%local_block%start(3):self%local_block%end(3) &
            )
        double precision, intent(inout) :: phi( &
            self%local_block%start(1):self%local_block%end(1), &
            self%local_block%start(2):self%local_block%end(2), &
            self%local_block%start(3):self%local_block%end(3) &
            )

        phi(:, :, :) = rho(:, :, :)

        block
            integer :: axis
            integer :: local_start(3), local_end(3)

            local_start = self%local_block%start
            local_end = self%local_block%end

            do axis = 1, 3
                if (self%local_block%start(axis) == self%global_block%start(axis)) then
                    phi(local_start(axis), :, :) = &
                        phi(local_start(axis), :, :) &
                        + self%boundary_condition_terms(1, axis)
                end if

                if (self%local_block%end(axis) == self%global_block%end(axis)) then
                    phi(local_end(axis), :, :) = &
                        phi(local_end(axis), :, :) &
                        + self%boundary_condition_terms(2, axis)
                end if
            end do
        end block

        block
            integer :: local_start(3), local_end(3)

            local_start = self%local_block%start
            local_end = self%local_block%end

            call self%fft3d%forward(phi(local_start(1):local_end(1), &
                                        local_start(2):local_end(2), &
                                        local_start(3):local_end(3)), &
                                    phi(local_start(1):local_end(1), &
                                        local_start(2):local_end(2), &
                                        local_start(3):local_end(3)))

            phi(local_start(1):local_end(1), &
                local_start(2):local_end(2), &
                local_start(3):local_end(3)) = &
                phi(local_start(1):local_end(1), &
                    local_start(2):local_end(2), &
                    local_start(3):local_end(3)) &
                /self%modified_wave_number(self%local_block%start(1):self%local_block%end(1), &
                                           self%local_block%start(2):self%local_block%end(2), &
                                           self%local_block%start(3):self%local_block%end(3))

            if (all(self%local_block%start == self%global_block%start)) then
                phi(local_start(1), &
                    local_start(2), &
                    local_start(3)) = 0d0
            end if

            call self%fft3d%backward(phi(local_start(1):local_end(1), &
                                         local_start(2):local_end(2), &
                                         local_start(3):local_end(3)), &
                                     phi(local_start(1):local_end(1), &
                                         local_start(2):local_end(2), &
                                         local_start(3):local_end(3)))
        end block
    end subroutine

end module
