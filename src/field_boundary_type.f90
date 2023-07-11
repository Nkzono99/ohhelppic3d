module m_field_boundary_type
    implicit none

    private
    public Field_BoundaryType_Periodic
    public Field_BoundaryType_Dirichlet
    public Field_BoundaryType_Neumann
    public Field_BoundaryType_Dirichlet_Neumann
    public Field_BoundaryType_Neumann_Dirichlet

    !> Periodic boundary type.
    !>
    !> Example (Array of length n = g.e-g.s, Logical size = n)
    !>  *----*----*--- ... ---*----*----o
    !> g.s  +1   +2          -2   -1   g.e(=g.s)
    !>
    !> g: global range (s: start, e: end)
    !>
    !> *: Actual data element (required)
    !> o: Actual data element (not required)
    integer, parameter :: Field_BoundaryType_Periodic = 0

    !> Dirichlet boundary type.
    !>
    !> Example (Array of length n = g.e-g.s-1, Logical size = 2(n+1))
    !>  D----*----*--- ... ---*----*----D
    !> g.s  +1   +2          -2   -1   g.e
    !>
    !> g: global range (s: start, e: end)
    !>
    !> *: Actual data element (required)
    !> D: Dirichlet boundary condition value (= 0)
    integer, parameter :: Field_BoundaryType_Dirichlet = 1

    !> Neumann boundary type.
    !>
    !> Example (Array of length n = g.e-g.s+1, Logical size = 2(n-1))
    !>  *N----*----*--- ... ---*-----*----*N
    !>  g.s  +1   +2          -2    -1   g.e
    !>
    !> g: global range (s: start, e: end)
    !>
    !> *: Actual data element (required)
    !> N: Neumann boundary condition value (= 0)
    integer, parameter :: Field_BoundaryType_Neumann = 2

    !> Dirichlet(left) and Neumann(right) boundary type. (Not tested to work properly.)
    !>
    !> Example (Array of length n = g.e-g.s, Logical size = 2n)
    !>  D----*----*--- ... ---*----*----*N
    !> g.s  +1   +2          -2   -1    g.e
    !>
    !> g: global range (s: start, e: end)
    !>
    !> *: Actual data element (required)
    !> D: Dirichlet boundary condition value (= 0)
    !> N: Neumann boundary condition value (= 0)
    integer, parameter :: Field_BoundaryType_Dirichlet_Neumann = 3

    !> Neumann(left) and Dirichlet(right) boundary type. (Not tested to work properly.)
    !>
    !> Example (Array of length n = g.e-g.s, Logical size = 2n)
    !>  *N----*----*--- ... ---*----*----D
    !>  g.s  +1   +2          -2   -1   g.e
    !>
    !> g: global range (s: start, e: end)
    !>
    !> *: Actual data element (required)
    !> D: Dirichlet boundary condition value (= 0)
    !> N: Neumann boundary condition value (= 0)
    integer, parameter :: Field_BoundaryType_Neumann_Dirichlet = 4

end module
