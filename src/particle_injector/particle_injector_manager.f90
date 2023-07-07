module m_particle_injector_manager
    use m_particle_injector
    use m_ohhelp
    use m_parameters
    use m_toml_wrapper, only: t_StringHolder
    use m_particle_injector_with_distribution, only: t_ParticleInjectorWithDistribution, &
                                                     new_ParticleInjectorWithDistribution
    use m_position_distribution
    use m_no_position_distribution
    use m_position_random_uniform_distribution
    use m_velocity_distribution
    use m_maxwellian_distribution
    use m_no_particle_injector

    use m_random_generator
    implicit none

    public t_ParticleInjectorManager, new_ParticleInjectorManager

    type t_ParticleInjectorHolder
        class(t_ParticleInjector), allocatable :: injector
    end type

    ! TODO: Change the name if come up with a name other than "Manager",
    type t_ParticleInjectorManager
        type(t_ParticleInjectorHolder), allocatable :: injectors_for_initialization(:)
        type(t_ParticleInjectorHolder), allocatable :: injectors_for_injection(:)
        class(t_RandomGenerator), pointer :: random_generator

    contains

        procedure :: inject_particles => particleInjectorManager_inject_particles
        procedure :: initialize_particles => particleInjectorManager_initialize_particles

    end type

contains

    function new_ParticleInjectorManager(parameters, random_generator) result(obj)
        type(t_ParticleInjectorManager) :: obj
        class(t_Parameters), intent(inout) :: parameters
        class(t_RandomGenerator), pointer, intent(in) :: random_generator

        integer :: ispec
        character(:), allocatable :: init_type

        obj%random_generator => random_generator

        allocate (obj%injectors_for_initialization(parameters%nspecies))

        do ispec = 1, parameters%nspecies
            init_type = parameters%plasma_initialization(ispec)%string

            select case (init_type)
            case ('random-uniform')
                block
                    class(t_PositionDistribution3d), allocatable :: pdist
                    class(t_VelocityDistribution3d), allocatable :: vdist

                    pdist = new_SimplePositionDistribution3d( &
                            new_PositionRandomUniformDistribution1d([0d0, 1d0*parameters%nx], obj%random_generator), &
                            new_PositionRandomUniformDistribution1d([0d0, 1d0*parameters%nx], obj%random_generator), &
                            new_PositionRandomUniformDistribution1d([0d0, 1d0*parameters%nx], obj%random_generator))

                    vdist = new_NoVelocityDistribution3d()
                    obj%injectors_for_initialization(ispec)%injector = &
                        new_ParticleInjectorWithDistribution(ispec, pdist, vdist, obj%random_generator)
                end block
            case ('none')
                block
                    obj%injectors_for_initialization(ispec)%injector = new_NoParticleInjector()
                end block
            end select
        end do

        ! TODO:実装
        block 
            integer :: nemission_types

            nemission_types = parameters%toml%require_int('plasma', 'injection', 'nemission_types')
            allocate (obj%injectors_for_injection(nemission_types))

        end block
    end function

    subroutine particleInjectorManager_initialize_particles(self, ohhelp)
        class(t_ParticleInjectorManager), intent(in) :: self
        class(t_OhHelp), intent(inout) :: ohhelp

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

        ! do i = 1, size(self%injectors_for_injection)
        !     call self%injectors_for_injection(i)%injector%inject_particles(1000, ohhelp)
        ! end do
    end subroutine

end module
