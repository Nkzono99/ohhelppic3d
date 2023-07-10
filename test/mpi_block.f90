module m_mpi_block
    use m_block
    implicit none

    private
    public create_local_block

contains

    subroutine create_local_block(global_block, nnodes, myid, local_block, ids)
        type(t_Block), intent(in) :: global_block
        integer, intent(in) :: nnodes(3)
        integer, intent(in) :: myid
        type(t_Block), intent(out) :: local_block
        integer, intent(out) :: ids(3)

        integer :: myid_x, myid_y, myid_z
        integer :: start(3), end(3)

        myid_x = mod(myid, nnodes(1))
        myid_y = mod(myid/nnodes(1), nnodes(2))
        myid_z = myid/nnodes(1)/nnodes(2)

        ids = [myid_x, myid_y, myid_z]

        start = ids*global_block%sizes/nnodes + global_block%start
        end = (ids + 1)*global_block%sizes/nnodes + global_block%start - 1
        if (myid_x == nnodes(1)) then
            end(1) = global_block%end(1)
        end if
        if (myid_y == nnodes(2)) then
            end(2) = global_block%end(2)
        end if
        if (myid_z == nnodes(3)) then
            end(3) = global_block%end(3)
        end if

        local_block = new_Block(start, end)
    end subroutine
end module
