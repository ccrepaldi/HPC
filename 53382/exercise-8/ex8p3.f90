 program bcast_ring
  use mpi
  implicit none
  integer,dimension(mpi_status_size) :: status
  integer, parameter :: tag=50
  integer :: id, ntasks, rc, i, msg, sum, dest_id, maxcount

  call mpi_init(rc)
  if (rc/=mpi_success) then
     print *,'MPI initialization failed.'
     stop
  end if
  call mpi_comm_size(MPI_COMM_WORLD,ntasks,rc)
  call mpi_comm_rank(MPI_COMM_WORLD,id,rc)

 ! print *, 'MPI OK, tasks: ', ntasks, 'ID: ', id

  if (mod(ntasks,2)/=0) then ! ntasks is odd
     if (id==0) print *, 'Processes must be a power of 2.'
     stop
  end if

  sum=id

  if (id/=0) then
     if (mod(id,2)/=0) then ! odd id
        dest_id=id-1
        call mpi_send(sum,1,mpi_integer,dest_id,tag,mpi_comm_world,rc)
        ! print *, "ID: ", id, "DEST: ", dest_id
     else ! even id
        if (id==ntasks-2) then
           maxcount=1
        else
           maxcount=2
        end if
        do i=1,maxcount
           call mpi_recv(msg,1,mpi_integer,mpi_any_source,tag,mpi_comm_world,status,rc)
           sum=sum+msg
        end do
        dest_id=id-2
        ! print *, "ID: ", id, "DEST: ", dest_id
        call mpi_send(sum,1,mpi_integer,dest_id,tag,mpi_comm_world,rc)
        ! print *, "ID: ", id, "SUM: ", sum
     end if
  else ! id is 0
     if (ntasks==2) then
        maxcount=1
     else
        maxcount=2
     end if
     do i=1,maxcount
        call mpi_recv(msg,1,mpi_integer,mpi_any_source,tag,mpi_comm_world,status,rc)
        sum=sum+msg
     end do
     print "(a,' ',i0,' is ',i0,'.')", trim('Sum from 0 to'), ntasks-1, sum
  end if
     
  call mpi_finalize(rc)

end program bcast_ring
