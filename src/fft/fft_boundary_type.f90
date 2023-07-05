module m_fft_boundary_type
    implicit none

    !> Periodic boundary type.
    !>
    !> Example (Array of length n, Logical size = n)
    !> *---*---*--- ... ---*---*---*---o
    !> 1   2   3          n-2 n-1  n  n+1(= 1)
    !>
    !> *: Actual data element (required)
    !> o: Actual data element (not required)
    integer, parameter :: FFT_BoundaryType_Periodic = 0

    !> Dirichlet boundary type.
    !> 
    !> Example (Array of length n, Logical size = 2(n+1))
    !> D---*---*---*--- ... ---*---*---*---D
    !> 0   1   2   3          n-2 n-1  n  n+1
    !>
    !> *: Actual data element (required)
    !> D: Dirichlet boundary condition value (= 0)
    integer, parameter :: FFT_BoundaryType_Dirichlet = 1

    !> Neumann boundary type.
    !>
    !> Example (Array of length n, Logical size = 2(n-1))
    !> *N---*---*--- ... ---*---*---*N
    !>  1   2   3          n-2 n-1  n
    !>
    !> *: Actual data element (required)
    !> N: Neumann boundary condition value (= 0)
    integer, parameter :: FFT_BoundaryType_Neumann = 2

    !> Dirichlet(left) and Neumann(right) boundary type. (Not tested to work properly.)
    !>
    !> Example (Array of length n, Logical size = 2n)
    !> D---*---*---*--- ... ---*---*---*N
    !> 0   1   2   3          n-2 n-1  n
    !>
    !> *: Actual data element (required)
    !> D: Dirichlet boundary condition value (= 0)
    !> N: Neumann boundary condition value (= 0)
    integer, parameter :: FFT_BoundaryType_Dirichlet_Neumann = 3

    !> Neumann(left) and Dirichlet(right) boundary type. (Not tested to work properly.)
    !>
    !> Example (Array of length n, Logical size = 2n)
    !> *N---*---*--- ... ---*---*---*---D
    !>  1   2   3          n-2 n-1  n  n+1
    !>
    !> *: Actual data element (required)
    !> D: Dirichlet boundary condition value (= 0)
    !> N: Neumann boundary condition value (= 0)
    integer, parameter :: FFT_BoundaryType_Neumann_Dirichlet = 4

    private
    public FFT_BoundaryType_Periodic
    public FFT_BoundaryType_Dirichlet
    public FFT_BoundaryType_Neumann
    public FFT_BoundaryType_Dirichlet_Neumann
    public FFT_BoundaryType_Neumann_Dirichlet

end module