module m_particle_boundary
    use oh_type
    implicit none

    private
    public t_ParticleBoundary

    type, abstract :: t_ParticleBoundary
    contains
        procedure(particleBoundary_apply), deferred :: apply
    end type

    interface
        subroutine particleBoundary_apply(self, particle, dt)
            import t_ParticleBoundary
            import oh_particle
            class(t_ParticleBoundary), intent(in) :: self
            type(oh_particle), intent(inout) :: particle
            double precision, intent(in) :: dt
        end subroutine
    end interface

contains

end module
