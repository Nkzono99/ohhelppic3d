module m_random_generator
    use m_science_constants, only: pi
    implicit none

    private
    public t_RandomGenerator

    !TODO: 乱数の取り扱いを実装
    type, abstract :: t_RandomGenerator
    contains

        procedure(randomGenerator_rand), deferred :: rand
        procedure :: normal => randomGenerator_normal
        procedure(randomGenerator_advance), deferred :: advance
        procedure :: random_fix => randomGenerator_random_fix
    end type

    interface
        function randomGenerator_rand(self) result(ret)
            import t_RandomGenerator
            class(t_RandomGenerator), intent(in) :: self
            double precision :: ret
        end function

        subroutine randomGenerator_advance(self, n)
            import t_RandomGenerator
            class(t_RandomGenerator), intent(in) :: self
            integer(8), intent(in) :: n
        end subroutine
    end interface

contains

    function randomGenerator_normal(self) result(ret)
        class(t_RandomGenerator), intent(in) :: self
        double precision :: ret

        double precision :: r1, r2

        r1 = self%rand()
        r2 = self%rand()
        ret = 2.0d0*sqrt(-2.0d0*log(r1))*cos(2.0d0*pi*r2)
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
