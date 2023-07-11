module m_const_velocity_distribution
    use m_velocity_distribution, only: t_VelocityDistribution1d
    implicit none

    private
    public t_ConstVelocityDistribution1d
    public new_ConstVelocityDistribution1d

    type, extends(t_VelocityDistribution1d) :: t_ConstVelocityDistribution1d
        double precision :: constant_value
    contains
        procedure :: sample => constVelocityDistribution1d_sample
    end type

contains

    function new_ConstVelocityDistribution1d(constant_value) result(obj)
        double precision, intent(in) :: constant_value
        type(t_ConstVelocityDistribution1d) :: obj

        obj%constant_value = constant_value
    end function

    function constVelocityDistribution1d_sample(self) result(ret)
        class(t_ConstVelocityDistribution1d), intent(in) :: self
        double precision :: ret

        ret = self%constant_value
    end function constVelocityDistribution1d_sample

end module
