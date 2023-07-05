module m_random
    implicit none

    !TODO: 乱数の取り扱いを実装

contains

    function random_uniform() result(ret)
        double precision :: ret

        ret = 0.5d0
    end function

    function random_normal() result(ret)
        double precision :: ret

        ret = 0.5d0
    end function

    function random_fix(val) result(ret)
        double precision, intent(in) :: val
        integer :: ret

        integer :: fixed

        fixed = int(val)

        if (random_uniform() <= fixed) then
            ret = fixed + 1
        else
            ret = fixed
        end if
    end function

end module
