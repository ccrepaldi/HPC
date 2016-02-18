program mpiexample
  use mpi
  implicit none
  !include 'mpif.h'
  integer,parameter :: tag = 50
  integer :: id,ntasks,source_id,dest_id,rc,i,nlen
  integer,dimension(mpi_status_size) :: status
  integer,dimension(30) :: msg
  character(len=MPI_MAX_PROCESSOR_NAME) :: host,arg1,arg2

  call mpi_init(rc)
  if (rc/=mpi_success) then
     print *,'MPI initialization failed.'
     stop
  end if

  call mpi_comm_size(mpi_comm_world,ntasks,rc)
  call mpi_comm_rank(mpi_comm_world,id,rc)
  call mpi_get_processor_name(host,nlen,rc)

  !call get_command_argument(1,arg1)
  !call get_command_argument(2,arg2)
  !if (id==0) print '(a,a,a,a)','Arguments : ',trim(arg1),',',trim(arg2)
  !print '(a,i0,a,i0)','id ',id,'/',ntasks
  
  msg(1)=id

  msg(2:30)=23

  dest_id=id+1
  source_id=id-1

  if (id.LT.ntasks-1) then
     call mpi_send(msg,30,mpi_integer,dest_id,tag,mpi_comm_world,rc)
     print *, "process id: ", id, "recipient id: ", dest_id
  end if

  if (id/=0) then
     call mpi_recv(msg,30,mpi_integer,source_id,tag,mpi_comm_world,status,rc)
     print *, "sender id: ", source_id, "msg(1): ", msg(1)
  end if

  !print '(/,a,a)','host: ',trim(host)
  call mpi_finalize(rc)
end program mpiexample
