program reduce
  use mpi
  implicit none
  integer, parameter :: ROOT=0
  integer :: id,ntasks,rc
  ! integer, dimension(mpi_status_size) :: status
  integer :: x,sum
  real(8) :: t0,t1,timeused

  call mpi_init(rc)
  if (rc/=mpi_success) then
     print *,'MPI initialization failed.'
     stop
  end if
  call mpi_comm_size(mpi_comm_world,ntasks,rc)
  call mpi_comm_rank(mpi_comm_world,id,rc)

   x=id

  t0=MPI_Wtime()
  call mpi_reduce(x,sum,1,mpi_integer,mpi_sum,0,mpi_comm_world,rc)
  t1=MPI_Wtime()
  timeused=t1-t0

  if (id==0) print *, ntasks, timeused
! if (id==0) print *, sum

  call mpi_finalize(rc)

end program reduce
