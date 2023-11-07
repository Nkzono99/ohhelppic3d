module m_flux_velocity_maxwell_distribution3d
    use m_velocity_distribution, only: t_VelocityDistribution3d
    use m_random_generator
    use m_science_constants, only: pi
    use m_vector
    implicit none

    private
    public t_FluxVelocityMaxwellDistribution3d
    public new_FluxVelocityMaxwellDistribution3d

    type, extends(t_VelocityDistribution3d) :: t_FluxVelocityMaxwellDistribution3d
        double precision :: velocity_variance
        double precision :: flow_direction(3)
        double precision :: flow_perp1(3)
        double precision :: flow_perp2(3)
        class(t_RandomGenerator), pointer :: random_generator
    contains
        procedure :: sample => fluxVelocityMaxwellDistribution3d_sample
    end type

contains

    function new_FluxVelocityMaxwellDistribution3d(velocity_variance, flow_direction, random_generator) result(obj)
        double precision, intent(in) :: velocity_variance
        double precision, intent(in) :: flow_direction(3)
        class(t_RandomGenerator), pointer :: random_generator
        type(t_FluxVelocityMaxwellDistribution3d) :: obj

        obj%velocity_variance = velocity_variance
        obj%random_generator => random_generator
        obj%flow_direction = normalized(flow_direction)

        obj%flow_perp1 = perp_vec(flow_direction)
        obj%flow_perp2 = cross(flow_direction, obj%flow_perp1)
    end function

    function fluxVelocityMaxwellDistribution3d_sample(self) result(ret)
        class(t_FluxVelocityMaxwellDistribution3d), intent(in) :: self
        double precision :: ret(3)

        double precision :: rands(3)

        double precision :: sigma
        double precision :: vperp(2)
        double precision :: vpara

        rands(1) = self%random_generator%rand()
        rands(2) = self%random_generator%rand()
        rands(3) = self%random_generator%rand()

        ! Assign maxwellian velocity.
        sigma = self%velocity_variance

        ! v in (-inf, +inf)
        ! Probability density function:
        !   f(v) = 1/sqrt(2*pi*sigma^2)*exp(-v^2/(2*sigma^2))
        vperp(1) = 2.0d0*sigma*sqrt(-2.0d0*log(rands(1)))*cos(2.0d0*PI*rands(2))
        vperp(2) = 2.0d0*sigma*sqrt(-2.0d0*log(rands(1)))*sin(2.0d0*PI*rands(2))

        ! v in (0, +inf)
        ! Probability density function:
        !   f(v) = 1/sigma^2*v*exp(-v^2/(2*sigma^2))
        !   F(v) = 1 - exp(-v^2/(2*sigma^2))
        ! v(3) = sqrt(-2*sigma*sigma*log(rands(3)))
        vpara = sqrt(-2*sigma*sigma*log(rands(3)))

        ret(:) = self%flow_direction*vpara &
                     + self%flow_perp1*vperp(1) &
                     + self%flow_perp2*vperp(2)
    end function

end module
