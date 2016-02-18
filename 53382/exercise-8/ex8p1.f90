program comline_bcast
  use mpi
  implicit none
  integer, parameter :: root_id=0
  integer :: id, ntasks, rc, i, narg, arg_size
  integer, allocatable :: msg(:)
  character(len=80) :: arg

  call mpi_init(rc)
  if (rc/=mpi_success) then
     print *,'MPI initialization failed.'
     stop
  end if
  call mpi_comm_size(MPI_COMM_WORLD,ntasks,rc)
  call mpi_comm_rank(MPI_COMM_WORLD,id,rc)

  narg=command_argument_count()
  allocate(msg(narg))

  if (id==root_id) then
     do i=1,narg
        call get_command_argument(i, arg)
        read(arg,*) msg(i)
     end do
  end if
  call mpi_bcast(narg,1,mpi_integer,root_id,mpi_comm_world,rc) ! broadcasts the number of arguments
  call mpi_bcast(msg,narg,mpi_integer,root_id,mpi_comm_world,rc) ! broadcasts the arguments
  do 
     if(id==i) then
        ! print "(50('-'))"
        print *, 'RANK: ',  id, msg(:)
        ! print "(50('-'))"
     end if
     i=i+1
     if (i==ntasks) then 
        exit
     end if
  end do
  
  call mpi_finalize(rc)

end program comline_bcast
