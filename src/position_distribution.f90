module m_position_distribution
    implicit none

    private
    public t_PositionDistribution1d
    public t_PositionDistribution3d
    public t_SimplePositionDistribution3d, new_SimplePositionDistribution3d

    type, abstract :: t_PositionDistribution1d
    contains
        procedure(distribution1d_sample), deferred :: sample
        procedure(distribution1d_subdomain_ratio), deferred :: subdomain_ratio
    end type

    interface
        function distribution1d_sample(self, subdomain_range) result(ret)
            import t_PositionDistribution1d
            class(t_PositionDistribution1d), intent(in) :: self
            double precision, intent(in) :: subdomain_range(2)
            double precision :: ret
        end function

        function distribution1d_subdomain_ratio(self, subdomain_range) result(ret)
            import t_PositionDistribution1d
            class(t_PositionDistribution1d), intent(in) :: self
            double precision, intent(in) :: subdomain_range(2)
            double precision :: ret
        end function
    end interface

    type, abstract :: t_PositionDistribution3d
    contains
        procedure(distribution3d_sample), deferred :: sample
        procedure(distribution3d_subdomain_ratio), deferred :: subdomain_ratio
    end type

    interface
        function distribution3d_sample(self, subdomain_range) result(ret)
            import t_PositionDistribution3d
            class(t_PositionDistribution3d), intent(in) :: self
            double precision, intent(in) :: subdomain_range(2, 3)
            double precision :: ret(3)
        end function

        function distribution3d_subdomain_ratio(self, subdomain_range) result(ret)
            import t_PositionDistribution3d
            class(t_PositionDistribution3d), intent(in) :: self
            double precision, intent(in) :: subdomain_range(2, 3)
            double precision :: ret
        end function
    end interface

    type, extends(t_PositionDistribution3d) :: t_SimplePositionDistribution3d
        class(t_PositionDistribution1d), allocatable :: distribution_x
        class(t_PositionDistribution1d), allocatable :: distribution_y
        class(t_PositionDistribution1d), allocatable :: distribution_z
    contains
        procedure :: sample => simplePositionDistribution3d_sample
        procedure :: subdomain_ratio => simplePositionDistribution3d_subdomain_ratio
    end type

contains

    function new_SimplePositionDistribution3d(distribution_x, distribution_y, distribution_z) result(obj)
        class(t_PositionDistribution1d), intent(in) :: distribution_x
        class(t_PositionDistribution1d), intent(in) :: distribution_y
        class(t_PositionDistribution1d), intent(in) :: distribution_z
        type(t_SimplePositionDistribution3d) :: obj

        obj%distribution_x = distribution_x
        obj%distribution_y = distribution_y
        obj%distribution_z = distribution_z
    end function

    function simplePositionDistribution3d_sample(self, subdomain_range) result(ret)
        class(t_SimplePositionDistribution3d), intent(in) :: self
        double precision, intent(in) :: subdomain_range(2, 3)
        double precision :: ret(3)

        ret(1) = self%distribution_x%sample(subdomain_range(:, 1))
        ret(2) = self%distribution_y%sample(subdomain_range(:, 2))
        ret(3) = self%distribution_z%sample(subdomain_range(:, 3))
    end function

    function simplePositionDistribution3d_subdomain_ratio(self, subdomain_range) result(ret)
        class(t_SimplePositionDistribution3d), intent(in) :: self
        double precision, intent(in) :: subdomain_range(2, 3)
        double precision :: ret

        ret = self%distribution_x%subdomain_ratio(subdomain_range(:, 1)) &
              *self%distribution_y%subdomain_ratio(subdomain_range(:, 2)) &
              *self%distribution_z%subdomain_ratio(subdomain_range(:, 3))
    end function

end module
