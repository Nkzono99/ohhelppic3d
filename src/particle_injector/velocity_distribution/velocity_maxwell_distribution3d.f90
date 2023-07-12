module m_velocity_maxwell_distribution3d
    use m_velocity_distribution, only: t_VelocityDistribution3d
    use m_random_generator
    use m_science_constants, only: pi
    use m_vector
    implicit none

    private
    public t_VelocityMaxwellDistribution3d
    public new_VelocityMaxwellDistribution3d

    type, extends(t_VelocityDistribution3d) :: t_VelocityMaxwellDistribution3d
        double precision :: velocity_variance
        double precision :: flow_direction(3)
        double precision :: flow_perp1(3)
        double precision :: flow_perp2(3)
        class(t_RandomGenerator), pointer :: random_generator
    contains
        procedure :: sample => VelocityMaxwellDistribution3d_sample
    end type

contains

    function new_VelocityMaxwellDistribution3d(velocity_variance, random_generator) result(obj)
        double precision, intent(in) :: velocity_variance
        class(t_RandomGenerator), pointer :: random_generator
        type(t_VelocityMaxwellDistribution3d) :: obj

        obj%velocity_variance = velocity_variance
        obj%random_generator => random_generator
    end function

    function VelocityMaxwellDistribution3d_sample(self) result(ret)
        class(t_VelocityMaxwellDistribution3d), intent(in) :: self
        double precision :: ret(3)

        double precision :: rands(4)

        double precision :: sigma

        rands(1) = self%random_generator%rand()
        rands(2) = self%random_generator%rand()
        rands(3) = self%random_generator%rand()
        rands(4) = self%random_generator%rand()

        ! Assign maxwellian velocity.
        sigma = self%velocity_variance

        ! v in (-inf, +inf)
        ! Probability density function:
        !   f(v) = 1/sqrt(2*pi*sigma^2)*exp(-v^2/(2*sigma^2))
        ret(1) = 2.0d0*sigma*sqrt(-2.0d0*log(rands(1)))*cos(2.0d0*PI*rands(2))
        ret(2) = 2.0d0*sigma*sqrt(-2.0d0*log(rands(1)))*sin(2.0d0*PI*rands(2))
        ret(3) = 2.0d0*sigma*sqrt(-2.0d0*log(rands(3)))*sin(2.0d0*PI*rands(4))

        ret = max(-sigma*6, min(ret, sigma*6))
        ! print *, ret
    end function

end module
