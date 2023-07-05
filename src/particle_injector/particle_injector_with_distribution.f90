module m_particle_injector_with_distribution
    use m_ohhelp, only: t_OhHelp, oh_particle
    use m_position_distribution, only: t_PositionDistribution3d
    use m_velocity_distribution, only: t_VelocityDistribution3d
    use m_random
    use m_particle_injector
    implicit none

    private
    public t_ParticleInjectorWithDistribution, new_ParticleInjectorWithDistribution

    type, extends(t_ParticleInjector) :: t_ParticleInjectorWithDistribution
        class(t_PositionDistribution3d), allocatable :: distribution_for_position
        class(t_VelocityDistribution3d), allocatable :: distribution_for_velocity

    contains
        procedure :: inject_particles => particleInjectorWithDistribution_inject_particles
    end type

contains

    function new_ParticleInjectorWithDistribution(ispec, &
                                                  distribution_for_position, &
                                                  distribution_for_velocity) result(obj)
        integer, intent(in) :: ispec
        class(t_PositionDistribution3d), intent(in) :: distribution_for_position
        class(t_VelocityDistribution3d), intent(in) :: distribution_for_velocity
        type(t_ParticleInjectorWithDistribution) :: obj

        obj%ispec = ispec
        obj%distribution_for_position = distribution_for_position
        obj%distribution_for_velocity = distribution_for_velocity
    end function

    subroutine particleInjectorWithDistribution_inject_particles(self, nparticles, ohhelp)
        class(t_ParticleInjectorWithDistribution), intent(in) :: self
        integer, intent(in) :: nparticles
        class(t_OhHelp), intent(inout) :: ohhelp

        call inject_particles(self%ispec, &
                              self%distribution_for_position, &
                              self%distribution_for_velocity, &
                              nparticles, &
                              ohhelp)
    end subroutine

    subroutine inject_particles(ispec, &
                                distribution_for_position, &
                                distribution_for_velocity, &
                                nparticles, &
                                ohhelp)
        integer, intent(in) :: ispec
        class(t_PositionDistribution3d), intent(in) :: distribution_for_position
        class(t_VelocityDistribution3d), intent(in) :: distribution_for_velocity
        integer, intent(in) :: nparticles
        class(t_OhHelp), intent(inout) :: ohhelp

        double precision :: subdomain_range(2, 3)
        double precision :: ratio
        integer :: iparticle

        print *, ohhelp%subdomain_id
        print *, ohhelp%subdomain_range
        print *, ohhelp%subdomain_id(1) + 1
        print *, ohhelp%subdomain_range(:, :, ohhelp%subdomain_id(1) + 1)
        subdomain_range(:, :) = 1d0*ohhelp%subdomain_range(:, :, ohhelp%subdomain_id(1) + 1)
        ratio = distribution_for_position%subdomain_ratio(subdomain_range)

        do iparticle = 1, random_fix(nparticles*ratio)
            block
                type(oh_particle) :: particle
                double precision :: position(3)
                double precision :: velocity(3)

                position(:) = distribution_for_position%sample(subdomain_range)
                velocity(:) = distribution_for_velocity%sample()
                particle%x = position(1)
                particle%y = position(2)
                particle%z = position(3)
                particle%vx = velocity(1)
                particle%vy = velocity(2)
                particle%vz = velocity(3)
                particle%nid = ohhelp%subdomain_id(1)
                particle%pid = 0
                particle%spec = ispec

                call ohhelp%inject_particle(particle)
            end block
        end do
    end subroutine

end module
