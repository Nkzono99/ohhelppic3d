
module m_pcg_generator
    use m_random_generator, only: t_RandomGenerator
    use m_pcg32
    implicit none

    private
    public t_PcgGenerator
    public new_PcgGenerator

    type, extends(t_RandomGenerator) :: t_PcgGenerator
        type(pcg_state_setseq_64) :: rng
        integer(8) :: seeds(2)
    contains
        procedure :: rand => pcgGenerator_rand
        procedure :: advance => pcgGenerator_advance
    end type

contains

    function new_PcgGenerator(seeds) result(obj)
        type(t_PcgGenerator) :: obj
        integer(8), intent(in) :: seeds(2)

        obj%seeds = seeds
        call pcg32_srandom_r(obj%rng, seeds(1), seeds(2))
    end function

    function pcgGenerator_rand(self) result(ret)
        class(t_PcgGenerator), intent(in) :: self
        double precision :: ret

        ret = pcg32_random_double_r(self%rng)
    end function

    subroutine pcgGenerator_advance(self, n)
        class(t_PcgGenerator), intent(in) :: self
        integer(8), intent(in) :: n

        call pcg32_advance_r(self%rng, n)
    end subroutine

end module
