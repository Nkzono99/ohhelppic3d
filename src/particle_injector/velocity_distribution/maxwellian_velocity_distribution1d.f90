module m_maxwellian_velocity_distribution1d
    use m_velocity_distribution, only: t_VelocityDistribution1d
    use m_random_generator
    implicit none

    private
    public t_MaxwellianVelocityDistribution1d
    public new_MaxwellianVelocityDistribution1d

    type, extends(t_VelocityDistribution1d) :: t_MaxwellianVelocityDistribution1d
        class(t_RandomGenerator), pointer :: random_generator
    contains
        procedure :: sample => maxwellianVelocityDistribution1d_sample
    end type

contains

    function new_MaxwellianVelocityDistribution1d(random_generator) result(obj)
        class(t_RandomGenerator), pointer :: random_generator
        type(t_MaxwellianVelocityDistribution1d) :: obj

        obj%random_generator => random_generator
    end function

    function maxwellianVelocityDistribution1d_sample(self) result(ret)
        class(t_MaxwellianVelocityDistribution1d), intent(in) :: self
        double precision :: ret

        ! Upper and lower of self%domains.
        double precision :: dl(3), du(3)
        double precision :: v(3)

        double precision :: sigma

        double precision :: PI

        integer :: icon
        double precision :: rands(6)

        PI = acos(-1.0d0)

        ! Assign maxwellian velocity.
        sigma = self%vel_variance

        ! v in (-inf, +inf)
        ! Probability density function:
        !   f(v) = 1/sqrt(2*pi*sigma^2)*exp(-v^2/(2*sigma^2))
        v(1) = 2.0d0*sigma*sqrt(-2.0d0*log(rands(4)))*cos(2.0d0*PI*rands(5))
        v(2) = 2.0d0*sigma*sqrt(-2.0d0*log(rands(4)))*sin(2.0d0*PI*rands(5))

        ! v in (0, +inf)
        ! Probability density function:
        !   f(v) = 1/sigma^2*v*exp(-v^2/(2*sigma^2))
        !   F(v) = 1 - exp(-v^2/(2*sigma^2))
        ! v(3) = sqrt(-2*sigma*sigma*log(rands(6)))
        v(3) = sqrt(-2*sigma*sigma*log(rands(6)))

        rands(4:6) = self%vel_direction*v(3) &
                     + self%vel_perp1*v(1) &
                     + self%vel_perp2*v(2)

        success = .true.
    end function
    

end module
