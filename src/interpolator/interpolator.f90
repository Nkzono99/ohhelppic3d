module m_interpolator
    use m_ohhelp, only: oh_particle
    use m_ohfield, only: t_OhField
    implicit none

    type, abstract :: t_Interpolator

    contains
        procedure(interpolator_interp), deferred :: interp

        procedure :: local_position_from(global_position)
    end type

    interface
        function interpolator_interp(self, particle, field) result(ret)
            import t_Interpolator
            import oh_particle
            import t_OhField
            class(t_Interpolator), intent(in) :: self
            class(oh_particle), intent(in) :: particle
            class(t_OhField), intent(in) :: field
            double precision :: ret(field%extension_info%nelements)
        end function
    end interface

contains

    function local_position_from(self, global_position, subdomain_range)
        class(t_Interpolator), intent(in) :: self
        double precision :: global_position(3)
        double precision :: subdomain_range(2, 3)

        
    end function

end module
