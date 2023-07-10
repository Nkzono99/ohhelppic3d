module m_periodic_particle_boundary
    use m_particle_boundary
    use oh_type, only: oh_particle
    implicit none

    private
    public new_PeriodicParticleBoundary

    type, extends(t_ParticleBoundary) :: t_PeriodicParticleBoundary
        integer :: axis
        integer :: n
    contains
        procedure :: apply => periodicParticleBoundary_apply
    end type

contains

    function new_PeriodicParticleBoundary(axis, n) result(obj)
        integer, intent(in) :: axis
        integer :: n
        type(t_PeriodicParticleBoundary) :: obj

        obj%axis = axis
        obj%n = n
    end function

    subroutine periodicParticleBoundary_apply(self, particle, dt)
        class(t_PeriodicParticleBoundary), intent(in) :: self
        type(oh_particle), intent(inout) :: particle
        double precision, intent(in) :: dt

        select case (self%axis)
        case (1)
            call apply_periodic(self%n, particle%x)
        case (2)
            call apply_periodic(self%n, particle%y)
        case (3)
            call apply_periodic(self%n, particle%z)
        end select
    end subroutine

    subroutine apply_periodic(n, p)
        integer, intent(in) :: n
        double precision, intent(inout) :: p

        if (p < 0) then
            p = n + p
        end if

        if (p >= n) then
            p = p - n
        end if
    end subroutine

end module
