! Try using nprocs=2
program mpi_testing
  use mpi
  implicit none
  integer,parameter :: root=0,tag=1,nmax=20 
  integer :: rc,myid,nprocs,req,&
       &status(mpi_status_size) 
  integer :: count,winner,i,index 
  integer,allocatable :: x(:) 
  integer :: ios,j,data(nmax)
  integer,parameter :: svalue=42 
  logical :: done

  call mpi_init(rc)
  call mpi_comm_rank(mpi_comm_world,myid,rc) 
  call mpi_comm_size(mpi_comm_world,nprocs,rc)
  count=0
  call mpi_irecv(winner,1,mpi_integer,&
       &mpi_any_source,tag,mpi_comm_world,req,rc)
 ! print *, myid, nprocs, nmax, nmax/nprocs
  allocate(x(nmax/nprocs))
  if (myid == root) then
    ! print *, 'i am root'
     ! gets data from file
    ! print *, 'opening'
     open(unit=1,file='prob4.data',iostat=ios,status='old')
     if (ios/=0) then
        write (6,'(a)')'*** Error in opening file prob4.data'
        stop 
     end if
     do j=1,nmax
        read(1,*,iostat=ios) data(j)
        if (ios<0) exit
     end do
     close(1)
    ! print *, 'closing'
     ! send pieces of data to all processes
  end if
  ! send pieces of data to all processes                                     
  call mpi_scatter(data,nmax/nprocs,mpi_integer,x,nmax/nprocs,mpi_integer,&
       &root, mpi_comm_world, rc)

  searchloop: do index=1,(nmax/nprocs)
     if (x(index)==svalue) then
        print *,myid,': found ',x(index), ', index:', index+(nmax/nprocs)*myid  
        do i=0,nprocs-1
           call mpi_send(myid,1,mpi_integer,i,tag,mpi_comm_world,rc)
        end do
        exit searchloop
     end if
     call mpi_test(req,done,status,rc) 
     if (done) then
        print *,myid,': lost to ',winner, ', index:', index+(nmax/nprocs)*myid
        exit searchloop 
     end if
     if (index==nmax/nprocs) print *, myid, ': reached end of the piece of data, index: ', index+(nmax/nprocs)*myid
  end do searchloop
  call mpi_finalize(rc) 
end  program mpi_testing
