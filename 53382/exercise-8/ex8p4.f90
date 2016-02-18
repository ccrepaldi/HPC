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

  if (id/=0) then
     msg=id*10
     do i=1,100
        call mpi_send(msg,1,mpi_integer,0,tag,mpi_comm_world,rc)
     end do
  else
     do i=1,(ntasks-1)*100
        call mpi_recv(msg,1,mpi_integer,mpi_any_source,tag,mpi_comm_world,status,rc)
        print "(a, i0)", trim('MESSAGE: '), msg
     end do
  end if

  call mpi_finalize(rc)

end program bcast_ring
