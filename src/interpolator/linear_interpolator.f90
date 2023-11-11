module m_linear_interpolator
    use m_interpolator
    use oh_type, only: oh_particle
    use m_ohfield, only: t_OhField
    implicit none

    private
    public t_LinearInterpolator
    public new_LinearInterpolator

    type, extends(t_Interpolator) :: t_LinearInterpolator
    contains
        procedure :: interp => linearInterpolator_interp
    end type

contains

    function new_LinearInterpolator() result(obj)
        type(t_LinearInterpolator) :: obj
    end function

    function linearInterpolator_interp(self, particle, field, ps) result(ret)
        class(t_LinearInterpolator), intent(in) :: self
        type(oh_particle), intent(in) :: particle
        class(t_OhField), intent(in) :: field
        integer, intent(in) :: ps
        double precision :: ret(field%nelements)

        double precision :: local_position(3)
        integer :: ilocal_position(3)
        double precision :: v(field%nelements, 2, 2, 2)

        local_position = field%to_local_position([particle%x, particle%y, particle%z], ps)
        ilocal_position = int(local_position)

        block
            integer :: i, j, k
            integer :: ix, iy, iz
            double precision :: val(field%nelements)

            do concurrent(i=1:2, j=1:2, k=1:2)
                ix = ilocal_position(1) + i - 1
                iy = ilocal_position(2) + j - 1
                iz = ilocal_position(3) + k - 1
                val = field%values(:, ix, iy, iz, ps)
                v(:, i, j, k) = val
            end do
        end block

        block
            integer :: i, j, k
            double precision :: r(3)
            double precision :: rxyz(2, 3)

            r = local_position - ilocal_position

            rxyz(1, :) = r(:)
            rxyz(2, :) = 1d0 - r(:)

            ret(:) = 0d0

            do i = 1, 2
            do j = 1, 2
            do k = 1, 2
                ret(:) = ret(:) + v(:, i, j, k)*rxyz(i, 1)*rxyz(j, 2)*rxyz(k, 3)
            end do
            end do
            end do
        end block
    end function

end module
