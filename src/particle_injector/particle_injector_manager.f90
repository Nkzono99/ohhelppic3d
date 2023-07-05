module m_particle_injector_manager
    use m_particle_injector
    use m_ohhelp
    use m_parameters
    implicit none

    public t_ParticleInjectorManager, new_ParticleInjectorManager

    type t_ParticleInjectorHolder
        class(t_ParticleInjector), allocatable :: injector
    end type

    ! TODO: Change the name if come up with a name other than "Manager", 
    type t_ParticleInjectorManager
        type(t_ParticleInjectorHolder), allocatable :: injectors_for_initialization(:)
        type(t_ParticleInjectorHolder), allocatable :: injectors_for_injection(:)

    contains

        procedure :: inject_particles => particleInjectorManager_inject_particles
        procedure :: initialize_particles => particleInjectorManager_initialize_particles

    end type

contains

    function new_ParticleInjectorManager() result(obj)
        type(t_ParticleInjectorManager) :: obj

        ! TODO: 長くなってもいいのでインジェクターの生成を行う
    end function

    subroutine particleInjectorManager_initialize_particles(self, ohhelp, parameters)
        class(t_ParticleInjectorManager), intent(in) :: self
        class(t_OhHelp), intent(inout) :: ohhelp
        class(t_Parameters), intent(in) :: parameters

        integer :: i

        do i = 1, size(self%injectors_for_initialization)
            call self%injectors_for_initialization(i)%injector%inject_particles(1000, ohhelp)
        end do
    end subroutine

    subroutine particleInjectorManager_inject_particles(self, dt, ohhelp)
        class(t_ParticleInjectorManager), intent(in) :: self
        double precision, intent(in) :: dt
        class(t_OhHelp), intent(inout) :: ohhelp

        integer :: i

        do i = 1, size(self%injectors_for_injection)
            call self%injectors_for_injection(i)%injector%inject_particles(1000, ohhelp)
        end do
    end subroutine

end module
