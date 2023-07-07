module m_random_generator
    implicit none

    !TODO: 乱数の取り扱いを実装
    type, abstract :: t_RandomGenerator
    contains

        procedure(randomGenerator_rand), deferred :: rand
        procedure :: normal => randomGenerator_normal
        procedure :: random_fix => randomGenerator_random_fix
    end type

    interface
        function randomGenerator_rand(self) result(ret)
            import t_RandomGenerator
            class(t_RandomGenerator), intent(in) :: self
            double precision :: ret
        end function
    end interface

contains

    function randomGenerator_normal(self) result(ret)
        class(t_RandomGenerator), intent(in) :: self
        double precision :: ret

        ! TODO: 実装
        ret = 0
    end function

    function randomGenerator_random_fix(self, value) result(ret)
        class(t_RandomGenerator), intent(in) :: self
        double precision, intent(in) :: value
        double precision :: ret

        integer :: fixed
        fixed = int(value)

        if (self%rand() <= value - fixed) then
            ret = fixed + 1
        else
            ret = fixed
        end if
    end function

end module
