module m_hdf5_for_ohfield
    use m_hdf5
    use m_ohfield
    use m_domain
    implicit none

    private
    public t_Hdf5ForOhfield
    public new_Hdf5ForOhfield

    type :: t_Hdf5ForOhfield
        type(t_HDF5File) :: hdf5
        type(t_HDF5Group) :: group

        character(:), allocatable :: filename
        character(:), allocatable :: group_name
        integer(kind=8) :: global_shape(3)
        integer(kind=8) :: global_offset(3)
    contains
        procedure :: write => hdf5ForOhfield_write
        procedure :: close => hdf5ForOhfield_close
    end type

contains

    function new_Hdf5ForOhfield(filename, group_name, global_shape, global_offset, comm) result(obj)
        character(*), intent(in) :: filename
        character(*), intent(in) :: group_name
        integer(kind=8), intent(in) :: global_shape(3)
        integer(kind=8), intent(in) :: global_offset(3)
        integer, intent(in) :: comm
        type(t_Hdf5ForOhfield) :: obj

        obj%filename = filename
        obj%group_name = group_name
        obj%global_shape = global_shape
        obj%global_offset = global_offset

        obj%hdf5 = new_HDF5File(filename, 'w', comm)
        obj%group = obj%hdf5%create_group(group_name)
    end function

    subroutine hdf5ForOhfield_write(self, dataset_name, ohfield, ps)
        class(t_Hdf5ForOhfield), intent(in) :: self
        character(*), intent(in) :: dataset_name
        class(t_OhField), intent(in) :: ohfield
        integer, intent(in) :: ps

        integer :: xl, xu, yl, yu, zl, zu
        integer :: status
        type(t_SubDomain3d) :: domain_info

        xl = ohfield%subdomain_range(1, 1, ps)
        xu = ohfield%subdomain_range(2, 1, ps)
        yl = ohfield%subdomain_range(1, 2, ps)
        yu = ohfield%subdomain_range(2, 2, ps)
        zl = ohfield%subdomain_range(1, 3, ps)
        zu = ohfield%subdomain_range(2, 3, ps)

        block
            integer(kind=8) :: local_shape(3)
            integer(kind=8) :: local_offset(3)
            local_shape(:) = ohfield%local_shape(ps)
            local_offset(:) = [xl, yl, zl]
            domain_info = new_SubDomain3d( &
                          local_shape, &
                          local_offset, &
                          self%global_shape, &
                          self%global_offset)
        end block

        call self%group%write_dataset(dataset_name, &
                                      ohfield%values(1, 0:xu-xl, 0:yu-yl, 0:zu-zl, ps), &
                                      domain_info, &
                                      status)
    end subroutine

    subroutine hdf5ForOhfield_close(self)
        class(t_Hdf5ForOhfield), intent(in) :: self

        integer :: status

        call self%group%close(status)
        call self%hdf5%close(status)
    end subroutine

end module
