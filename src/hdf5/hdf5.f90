module m_hdf5
    use mpi
    use HDF5
    use m_domain

    implicit none

    private
    public hdf5_initialize, hdf5_finalize
    public t_HDF5Group
    public t_HDF5File, new_HDF5File

    type :: t_HDF5Group
        character(:), allocatable :: name !! Group name
        integer(kind=HID_T) :: id !! Group id
        character :: mode !! Access mode (r: readonly, w: writeonly, a: read/write)
    contains
        procedure, private :: write_double3d => hdf5Group_write_dataset3d_double
        generic :: write_dataset => write_double3d
        procedure :: create_group => hdf5Group_create_group
        procedure :: close => hdf5Group_close
    end type

    type, extends(t_HDF5Group) :: t_HDF5File
        character(:), allocatable :: filename !! File name

    contains
        procedure :: close => hdf5File_close
    end type

contains

    subroutine hdf5_initialize(status)
        integer, intent(out) :: status

        call h5open_f(status)
        call h5eset_auto_f(0, status)
    end subroutine

    subroutine hdf5_finalize(status)
        integer, intent(out) :: status

        call h5close_f(status)
    end subroutine

    function new_HDF5File(filename, mode, comm) result(obj)
        character(*), intent(in) :: filename !! File name to create or open
        character, intent(in) :: mode !! Access mode (r: readonly, w: writeonly, a: read/write)
        integer, intent(in) :: comm !! MPI communicator
        type(t_HDF5File) :: obj

        integer :: hdferr !! Error status
        integer(kind=HID_T):: prp_id !! Property list id

        obj%mode = mode
        obj%filename = filename
        obj%name = filename

        ! Configure property list for MPI parallel access.
        call h5pcreate_f(H5P_FILE_ACCESS_F, prp_id, hdferr)
        call h5pset_fapl_mpio_f(prp_id, comm, MPI_INFO_NULL, hdferr)

        if (mode == 'r') then
            call h5fopen_f(filename, H5F_ACC_RDONLY_F, obj%id, hdferr, access_prp=prp_id)
        else if (mode == 'w') then
            call h5fcreate_f(filename, H5F_ACC_TRUNC_F, obj%id, hdferr, access_prp=prp_id)
        else if (mode == 'a') then
            call h5fopen_f(filename, H5F_ACC_RDWR_F, obj%id, hdferr, access_prp=prp_id)
        end if

        call h5pclose_f(prp_id, hdferr)
    end function

    function hdf5Group_create_group(self, name) result(group)
        class(t_HDF5Group), intent(in) :: self
        character(*), intent(in) :: name
        type(t_HDF5Group) :: group

        integer :: hdferr !! Error status

        group%name = name
        group%mode = self%mode

        if (self%mode == 'w') then
            call h5gcreate_f(self%id, name, group%id, hdferr)
        else
            call h5gopen_f(self%id, name, group%id, hdferr)
        end if
    end function

    subroutine hdf5Group_read_dataset1d_int(self, dataset_name, dataset, status)
        class(t_HDF5Group), intent(in) :: self
        character(*), intent(in) :: dataset_name
        integer, intent(out) :: dataset(:)
        integer, intent(out) :: status

        integer, parameter :: ndims = 1
        integer(kind=HSIZE_T) :: dims(ndims)

        integer :: dataset_ndim
        integer(kind=HSIZE_T) :: dataset_dims(ndims), dataset_maxdims(ndims)

        integer(kind=HID_T) :: dataset_id
        integer(kind=HID_T) :: space_id

        integer :: idim

        dims = shape(dataset)

        call h5dopen_f(self%id, dataset_name, dataset_id, status)
        if (status /= 0) then
            print *, 'InvalidDatasetNameError: dataset name is invalid:', dataset_name
            return
        end if

        call h5dget_space_f(dataset_id, space_id, status)
        call h5sget_simple_extent_ndims_f(space_id, dataset_ndim, status)
        if (dataset_ndim /= ndims) then
            print *, 'InvalidDatasetRankError: dataset rank is invalid:', dataset_ndim
            return
        end if

        call h5sget_simple_extent_dims_f(space_id, dataset_dims, dataset_maxdims, status)
        do idim = 1, ndims
            if (dataset_dims(idim) > dims(idim)) then
                print *, 'InvalidDatasetDimensionError: dataset dimension is too large:', dataset_dims
                return
            end if
        end do

        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, dataset, dims, status)
        if (status /= 0) then
            return
        end if

        call h5sclose_f(space_id, status)
        call h5dclose_f(dataset_id, status)
    end subroutine

    subroutine hdf5Group_write_dataset3d_double(self, dataset_name, dataset, subdomain_info, status)
        class(t_HDF5Group), intent(in) :: self
        character(*), intent(in) :: dataset_name
        double precision, intent(in) :: dataset(:, :, :)
        type(t_SubDomain3d), intent(in) :: subdomain_info
        integer, intent(out) :: status

        integer, parameter :: ndims = 3

        integer(kind=HID_T) :: memory_space_id
        integer(kind=HID_T) :: file_space_id
        integer(kind=HID_T) :: dataset_id
        integer(kind=HID_T) :: property_id

        ! Create global dataset
        call h5screate_simple_f(ndims, subdomain_info%global_shape, file_space_id, status)
        call h5dcreate_f(self%id, dataset_name, H5T_NATIVE_DOUBLE, file_space_id, dataset_id, status)
        call h5sclose_f(file_space_id, status)

        ! Create local memory_space
        call h5screate_simple_f(ndims, subdomain_info%local_shape, memory_space_id, status)
        ! Create global file_space
        call h5dget_space_f(dataset_id, file_space_id, status)
        ! Convert local file_space
        call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, subdomain_info%local_offset, subdomain_info%local_shape, status)

        ! Set MPI collective
        call h5pcreate_f(H5P_DATASET_XFER_F, property_id, status)
        call h5pset_dxpl_mpio_f(property_id, H5FD_MPIO_COLLECTIVE_F, status)

        call h5dwrite_f(dataset_id, &
                        H5T_NATIVE_DOUBLE, &
                        dataset(1:subdomain_info%local_shape(1), 1:subdomain_info%local_shape(2), 1:subdomain_info%local_shape(3)), &
                        subdomain_info%global_shape, &
                        status, &
                        file_space_id=file_space_id, &
                        mem_space_id=memory_space_id, &
                        xfer_prp=property_id)

        call h5pclose_f(property_id, status)

        call h5sclose_f(file_space_id, status)
        call h5sclose_f(memory_space_id, status)

        call h5dclose_f(dataset_id, status)
    end subroutine

    subroutine hdf5Group_close(self, status)
        class(t_HDF5Group), intent(in) :: self
        integer, intent(out) :: status

        call h5gclose_f(self%id, status)
    end subroutine

    subroutine hdf5File_close(self, status)
        class(t_HDF5File), intent(in) :: self
        integer, intent(out) :: status

        call h5fclose_f(self%id, status)
    end subroutine

end module
