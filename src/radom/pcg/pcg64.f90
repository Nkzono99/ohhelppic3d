module m_pcg32
    use iso_c_binding, only: c_int64_t, c_double
    implicit none

    type, bind(c) :: pcg_state_setseq_64
        integer(c_int64_t) state; 
        integer(c_int64_t) inc; 
    end type

    interface
        subroutine pcg32_srandom_r(rng, state, initseq) bind(c, name='pcg_setseq_64_srandom_r')
            import pcg_state_setseq_64
            import c_int64_t
            type(pcg_state_setseq_64) :: rng
            integer(c_int64_t), value, intent(in) :: state
            integer(c_int64_t), value, intent(in) :: initseq
        end subroutine

        subroutine pcg32_advance_r(rng, delta) bind(c, name='pcg_setseq_64_advance_r')
            import pcg_state_setseq_64
            import c_int64_t
            type(pcg_state_setseq_64) :: rng
            integer(c_int64_t), value, intent(in) :: delta
        end subroutine

        function pcg32_random_double_r(rng) result(ret) bind(c)
            import pcg_state_setseq_64
            import c_double
            type(pcg_state_setseq_64) :: rng
            real(c_double) :: ret
        end function
    end interface

end module
