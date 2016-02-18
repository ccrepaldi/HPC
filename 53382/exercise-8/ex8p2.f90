program bcast_ring
  use mpi
  implicit none
  integer, parameter :: tag=50
  integer,dimension(mpi_status_size) :: status
  integer :: id, ntasks, rc, i, msg, dest_id, source_id

  call mpi_init(rc)
  if (rc/=mpi_success) then
     print *,'MPI initialization failed.'
     stop
  end if
  call mpi_comm_size(MPI_COMM_WORLD,ntasks,rc)
  call mpi_comm_rank(MPI_COMM_WORLD,id,rc)

  if (id==0) then
     msg=42
     call mpi_send(msg,1,mpi_integer,1,tag,mpi_comm_world,rc)
     print *, "ID: ", id, "MSG: ", msg
  else if (id.LT.ntasks-1) then
     source_id=id-1
     dest_id=id+1
     call mpi_recv(msg,1,mpi_integer,source_id,tag,mpi_comm_world,status,rc)
     call mpi_send(msg,1,mpi_integer,dest_id,tag,mpi_comm_world,rc)
     print *, "ID: ", id, "MSG: ", msg
  else
     source_id=id-1
     call mpi_recv(msg,1,mpi_integer,source_id,tag,mpi_comm_world,status,rc)
    print *, "ID: ", id, "MSG: ", msg
  end if
  call mpi_finalize(rc)

end program bcast_ring
