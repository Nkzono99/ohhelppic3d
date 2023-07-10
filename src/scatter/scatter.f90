module m_scatter
    use oh_type, only: oh_particle
    use m_ohfield, only: t_OhField
    implicit none

    private
    public t_Scatter

    type, abstract :: t_Scatter
    contains
        procedure(scatter_scatter), deferred :: scatter
    end type

    interface
        subroutine scatter_scatter(self, particle, ohfield, amount, ps)
            import t_Scatter
            import oh_particle
            import t_OhField
            class(t_Scatter), intent(in) :: self
            type(oh_particle), intent(in) :: particle
            class(t_OhField), intent(inout) :: ohfield
            double precision, intent(in) :: amount(ohfield%nelements)
            integer, intent(in) :: ps
        end subroutine
    end interface

contains

end module
