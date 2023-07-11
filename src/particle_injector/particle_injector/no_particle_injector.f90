module m_no_particle_injector
    use m_ohhelp
    use m_particle_injector
    implicit none

    private
    public t_NoParticleInjector
    public new_NoParticleInjector

    type, extends(t_ParticleInjector) :: t_NoParticleInjector
    contains
        procedure :: inject_particles => noParticleInjector_inject_particles
    end type

contains

    function new_NoParticleInjector() result(obj)
        type(t_NoParticleInjector) :: obj
    end function

    subroutine noParticleInjector_inject_particles(self, nparticles, ohhelp)
        class(t_NoParticleInjector), intent(in) :: self
        integer, intent(in) :: nparticles
        class(t_OhHelp), intent(inout) :: ohhelp
    end subroutine
    
    subroutine particleInjector_inject_particles(self, nparticles, ohhelp)
        class(t_ParticleInjector), intent(in) :: self
        integer, intent(in) :: nparticles
        class(t_OhHelp), intent(inout) :: ohhelp
    end subroutine

end module
