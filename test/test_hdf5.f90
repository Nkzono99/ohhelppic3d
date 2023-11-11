program test_hdf5
    use mpi
    use m_hdf5
    use m_domain
    implicit none

    character(*), parameter :: filename = 'example2.h5'
    integer :: ierr, myid, nprocs

    type(t_HDF5File) :: h5file
    type(t_HDF5Group) :: group
    double precision :: data3d(3, 2, 1)
    integer :: status
    type(t_SubDomain3d) :: domain_info

    call mpi_init(ierr)
    if (ierr /= 0) error stop "mpi_init failed"

    call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)

    if (ierr /= 0) error stop "mpi_comm_rank failed"

    call hdf5_initialize(status)

    data3d(1, 1, 1) = 1d-4*myid

    block
        integer(kind=8) :: local_shape(3)
        integer(kind=8) :: local_offset(3)
        integer(kind=8) :: global_shape(3)
        integer(kind=8) :: global_offset(3)
        local_shape(:) = [3, 2, 1]
        if (myid == 1) then
            local_offset(:) = [myid*3+1, 0, 0]
        else
            local_offset(:) = [myid*3, 0, 0]
        end if
        global_shape(:) = [nprocs*3, 2, 1]
        global_offset(:) = [0, 0, 0]
        domain_info = new_SubDomain3d( &
                      local_shape, &
                      local_offset, &
                      global_shape, &
                      global_offset &
                      )
    end block

    h5file = new_HDF5File(filename, 'w', MPI_COMM_WORLD)
    group = h5file%create_group('group1')
    call group%write_dataset('dataset1', data3d, domain_info, status)

    call group%close(status)
    call h5file%close(status)

    call hdf5_finalize(status)

    call mpi_finalize(ierr)
end program test_hdf5
