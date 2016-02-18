program always_works
  use mpi
  implicit none
  integer, parameter :: ROOT=0,tag=50
  integer :: id,ntasks,ierr,sbuf,rbuf
  integer, dimension(mpi_status_size) :: status

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,ntasks,ierr)
  call mpi_comm_rank(mpi_comm_world,id,ierr)

  sbuf=id
  
  ! always works
  if (id==0) then
     call mpi_send(sbuf,1,mpi_integer,1,tag,mpi_comm_world,ierr)
     call mpi_recv(rbuf,1,mpi_integer,1,mpi_any_tag,mpi_comm_world,status,ierr)
  else if (id==1) then
     call mpi_recv(rbuf,1,mpi_integer,0,mpi_any_tag,mpi_comm_world,status,ierr)
     call mpi_send(sbuf,1,mpi_integer,0,tag,mpi_comm_world,ierr)
  end if

  print *, ntasks, id, rbuf

  call mpi_finalize(ierr)

end program always_works
