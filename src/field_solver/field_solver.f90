module m_field_solver
    use m_ohhelp
    use m_ohfield
    implicit none

    private
    public t_FieldSolver

    type, abstract :: t_FieldSolver
    contains
        procedure(fieldSolver_solve), deferred :: solve
    end type

    interface
        subroutine fieldSolver_solve(self, dt, rho, aj, eb, phi, ohhelp)
            import t_FieldSolver
            import t_OhField
            import t_OhHelp
            class(t_FieldSolver), intent(in) :: self
            double precision, intent(in) :: dt
            class(t_OhField), intent(in) :: rho
            class(t_OhField), intent(in) :: aj
            class(t_OhField), intent(inout) :: eb
            class(t_OhField), intent(inout) :: phi
            class(t_OhHelp), intent(inout) :: ohhelp
        end subroutine
    end interface

contains

end module
