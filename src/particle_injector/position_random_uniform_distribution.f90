module m_position_random_uniform_distribution
    use m_position_distribution, only: t_PositionDistribution1d
    use m_random
    implicit none

    private
    public t_PositionRandomUniformDistribution1d
    public new_PositionRandomUniformDistribution1d

    type, extends(t_PositionDistribution1d) :: t_PositionRandomUniformDistribution1d
        double precision :: range(2)
        double precision, private :: total_range
    contains
        procedure :: sample => positionRandomUniformDistribution1d_sample
        procedure :: subdomain_ratio => positionRandomUniformDistribution1d_subdomain_range
    end type

contains

    function new_PositionRandomUniformDistribution1d(range) result(obj)
        type(t_PositionRandomUniformDistribution1d) :: obj
        double precision, intent(in) :: range(2)

        obj%range(:) = range(:)
        obj%total_range = range(2) - range(1)
    end function

    function positionRandomUniformDistribution1d_sample(self, subdomain_range) result(ret)
        class(t_PositionRandomUniformDistribution1d), intent(in) :: self
        double precision, intent(in) :: subdomain_range(2)
        double precision :: ret

        double precision :: lower, upper

        upper = min(self%range(2), subdomain_range(2))
        lower = max(self%range(1), subdomain_range(1))

        ret = random_uniform()*(upper - lower) + lower
    end function

    function positionRandomUniformDistribution1d_subdomain_range(self, subdomain_range) result(ret)
        class(t_PositionRandomUniformDistribution1d), intent(in) :: self
        double precision, intent(in) :: subdomain_range(2)
        double precision :: ret

        double precision :: lower, upper

        upper = min(self%range(2), subdomain_range(2))
        lower = max(self%range(1), subdomain_range(1))

        ret = max(upper - lower, 0d0)/(self%total_range)
    end function

end module