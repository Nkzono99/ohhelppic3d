module m_velocity_distribution
    type, abstract :: t_VelocityDistribution1d
    contains
        procedure(distribution1d_sample), deferred :: sample
    end type

    interface
        function distribution1d_sample(self) result(ret)
            import t_VelocityDistribution1d
            class(t_VelocityDistribution1d), intent(in) :: self
            double precision :: ret
        end function
    end interface

    type, abstract :: t_VelocityDistribution3d
    contains
        procedure(distribution3d_sample), deferred :: sample
    end type

    interface
        function distribution3d_sample(self) result(ret)
            import t_VelocityDistribution3d
            class(t_VelocityDistribution3d), intent(in) :: self
            double precision :: ret(3)
        end function
    end interface

    type, extends(t_VelocityDistribution3d) :: t_NoVelocityDistribution3d
    contains
        procedure :: sample => noVelocityDistribution3d_sample
    end type

    type, extends(t_VelocityDistribution3d) :: t_SimpleVelocityDistribution3d
        class(t_VelocityDistribution1d), allocatable :: distribution_x
        class(t_VelocityDistribution1d), allocatable :: distribution_y
        class(t_VelocityDistribution1d), allocatable :: distribution_z
    contains
        procedure :: sample => simpleVelocityDistribution3d_sample
    end type

contains

    function new_NoVelocityDistribution3d() result(obj)
        type(t_NoVelocityDistribution3d) :: obj

    end function

    function noVelocityDistribution3d_sample(self) result(ret)
        class(t_NoVelocityDistribution3d), intent(in) :: self
        double precision :: ret(3)

        ret(:) = 1.0d0
    end function

    function new_SimpleVelocityDistribution3d(distribution_x, distribution_y, distribution_z) result(obj)
        class(t_VelocityDistribution1d), intent(in) :: distribution_x
        class(t_VelocityDistribution1d), intent(in) :: distribution_y
        class(t_VelocityDistribution1d), intent(in) :: distribution_z
        type(t_SimpleVelocityDistribution3d) :: obj

        obj%distribution_x = distribution_x
        obj%distribution_y = distribution_y
        obj%distribution_z = distribution_z
    end function

    function simpleVelocityDistribution3d_sample(self) result(ret)
        class(t_SimpleVelocityDistribution3d), intent(in) :: self
        double precision :: ret(3)

        ret(1) = self%distribution_x%sample()
        ret(2) = self%distribution_y%sample()
        ret(3) = self%distribution_z%sample()
    end function

end module
