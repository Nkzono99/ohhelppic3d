module m_linear_scatter
    use m_scatter
    use oh_type, only: oh_particle
    use m_ohfield, only: t_OhField
    implicit none

    private
    public t_LinearScatter
    public new_LinearScatter

    type, extends(t_Scatter) :: t_LinearScatter
    contains
        procedure :: scatter => linearScatter_scatter
    end type

contains

    function new_LinearScatter() result(obj)
        type(t_LinearScatter) :: obj
    end function

    subroutine linearScatter_scatter(self, particle, ohfield, amount, ps)
        class(t_LinearScatter), intent(in) :: self
        type(oh_particle), intent(in) :: particle
        class(t_OhField), intent(inout) :: ohfield
        double precision, intent(in) :: amount(ohfield%nelements)
        integer, intent(in) :: ps

        double precision :: local_position(3)
        double precision :: ilocal_position(3)

        local_position = ohfield%to_local_position([particle%x, particle%y, particle%z], ps)
        ilocal_position = int(local_position)

        block
            integer :: i, j, k
            double precision :: r(3)
            double precision :: rxyz(2, 3)
            integer :: ix, iy, iz

            r = local_position - ilocal_position

            rxyz(1, :) = r(:)
            rxyz(2, :) = 1d0 - r(:)

            do concurrent(i=1:2, j=1:2, k=1:2)
                ix = ilocal_position(1) + i - 1
                iy = ilocal_position(2) + j - 1
                iz = ilocal_position(3) + k - 1
                ohfield%values(:, ix, iy, iz, ps) = &
                    ohfield%values(:, ix, iy, iz, ps) + &
                    amount*rxyz(i, 1)*rxyz(j, 2)*rxyz(k, 3)
            end do
        end block
    end subroutine

end module
