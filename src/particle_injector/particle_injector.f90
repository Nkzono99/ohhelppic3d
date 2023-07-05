module m_particle_injector
    use m_ohhelp, only: t_OhHelp
    use m_random
    implicit none

    private
    public t_ParticleInjector

    type, abstract :: t_ParticleInjector
        integer :: ispec
    contains
        procedure(particleInjector_inject_particles), deferred :: inject_particles
    end type

    interface
        subroutine particleInjector_inject_particles(self, nparticles, ohhelp)
            import t_ParticleInjector
            import t_OhHelp
            class(t_ParticleInjector), intent(in) :: self
            integer, intent(in) :: nparticles
            class(t_OhHelp), intent(inout) :: ohhelp
        end subroutine
    end interface

end module
