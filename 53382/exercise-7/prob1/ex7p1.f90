program mpiexample
  use mpi
  implicit none
  !include 'mpif.h'
  integer,parameter :: tag = 50
  integer :: id,ntasks,source_id,dest_id,rc,i,nlen
  integer,dimension(mpi_status_size) :: status
  integer,dimension(2) :: msg
  character(len=MPI_MAX_PROCESSOR_NAME) :: host,arg1,arg2

  call mpi_init(rc)
  if (rc/=mpi_success) then
     print *,'MPI initialization failed.'
     stop
  end if

  call mpi_comm_size(mpi_comm_world,ntasks,rc)
  call mpi_comm_rank(mpi_comm_world,id,rc)
  call mpi_get_processor_name(host,nlen,rc)

  call get_command_argument(1,arg1)
  call get_command_argument(2,arg2)
  if (id==0) print '(a,a,a,a)','Arguments : ',trim(arg1),',',trim(arg2)
  print '(a,i0,a,i0)','id ',id,'/',ntasks
  

  if (id/=0) then
     msg(1)=id
     msg(2)=ntasks
     dest_id=0
     call mpi_send(msg,2,mpi_integer,dest_id,tag,mpi_comm_world,rc)
  else
     do i=1,ntasks-1
        call mpi_recv(msg,2,mpi_integer,mpi_any_source,tag,mpi_comm_world,status,rc)
        source_id=status(mpi_source)
        print '(/,a,2(x,i0),5x,a,10(x,i0))','message:',msg,'sender:',status
     end do
  end if

  print '(/,a,a)','host: ',trim(host)
  call mpi_finalize(rc)
end program mpiexample
