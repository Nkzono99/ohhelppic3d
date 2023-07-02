module m_particle_mover
    use m_ohhelp, only: oh_particle
    implicit none

    private
    public t_ParticleMover

    type, abstract :: t_ParticleMover
    contains

        procedure(particleMover_move), deferred :: move

    end type

    interface
        subroutine particleMover_move(self, particle, qm, eb, dt)
            import t_ParticleMover
            import oh_particle
            class(t_ParticleMover), intent(in) :: self
            type(oh_particle), intent(inout) :: particle
            double precision, intent(in) :: qm
            double precision, intent(in) :: eb(6)
            double precision :: dt
        end subroutine
    end interface

end module
