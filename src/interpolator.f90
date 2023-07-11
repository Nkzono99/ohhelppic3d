module m_interpolator
    use oh_type, only: oh_particle
    use m_ohfield, only: t_OhField
    implicit none

    private
    public t_Interpolator

    type, abstract :: t_Interpolator

    contains
        procedure(interpolator_interp), deferred :: interp
    end type

    interface
        function interpolator_interp(self, particle, field, ps) result(ret)
            import t_Interpolator
            import oh_particle
            import t_OhField
            class(t_Interpolator), intent(in) :: self
            type(oh_particle), intent(in) :: particle
            class(t_OhField), intent(in) :: field
            integer, intent(in) :: ps
            double precision :: ret(field%nelements)
        end function
    end interface

end module
