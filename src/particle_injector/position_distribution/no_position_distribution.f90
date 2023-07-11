module m_no_position_distribution
    use m_position_distribution, only: t_PositionDistribution3d
    implicit none

    private
    public t_NoPositionDistribution3d, new_NoPositionDistribution3d

    type, extends(t_PositionDistribution3d) :: t_NoPositionDistribution3d
    contains
        procedure :: sample => noPositionDistribution3d_sample
        procedure :: subdomain_ratio => noPositionDistribution3d_subdomain_ratio
    end type

contains

    function new_NoPositionDistribution3d() result(obj)
        type(t_NoPositionDistribution3d) :: obj
    end function

    function noPositionDistribution3d_sample(self, subdomain_range) result(ret)
        class(t_NoPositionDistribution3d), intent(in) :: self
        double precision, intent(in) :: subdomain_range(2, 3)
        double precision :: ret(3)

        ret(:) = 0d0
    end function

    function noPositionDistribution3d_subdomain_ratio(self, subdomain_range) result(ret)
        class(t_NoPositionDistribution3d), intent(in) :: self
        double precision, intent(in) :: subdomain_range(2, 3)
        double precision :: ret

        ret = 0d0
    end function

end module
