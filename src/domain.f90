module m_domain

    type :: t_SubDomain3d
        integer(kind=8) :: local_shape(3)
        integer(kind=8) :: local_offset(3)
        integer(kind=8) :: global_shape(3)
        integer(kind=8) :: global_offset(3)
    end type

    private
    public t_SubDomain3d, new_SubDomain3d

contains

    function new_SubDomain3d(local_shape, local_offset, global_shape, global_offset) result(obj)
        integer(kind=8), intent(in) :: local_shape(3)
        integer(kind=8), intent(in) :: local_offset(3)
        integer(kind=8), intent(in) :: global_shape(3)
        integer(kind=8), intent(in) :: global_offset(3)

        type(t_SubDomain3d) :: obj

        obj%local_shape = local_shape
        obj%local_offset = local_offset
        obj%global_shape = global_shape
        obj%global_offset = global_offset
    end function

end module
