program probably_works
  use mpi
  implicit none
  integer, parameter :: ROOT=0,tag=50
  integer :: i,id,ntasks,ierr,msg_size
  integer*4, allocatable :: sbuf(:),rbuf(:) ! 32 bits integer arrays
  integer, dimension(mpi_status_size) :: status

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,ntasks,ierr)
  call mpi_comm_rank(mpi_comm_world,id,ierr)

  do i=1,100000 ! try with different message sizes
     if (allocated(rbuf)) then
        deallocate(rbuf,sbuf)
     end if

     allocate(sbuf(i),rbuf(i))

     ! probably works
     if (id==0) then
        call mpi_send(sbuf,i,mpi_integer,1,tag,mpi_comm_world,ierr)
        call mpi_recv(rbuf,i,mpi_integer,1,mpi_any_tag,mpi_comm_world,&
             &status,ierr)
     else if (id==1) then
        call mpi_send(sbuf,i,mpi_integer,0,tag,mpi_comm_world,ierr)
        call mpi_recv(rbuf,i,mpi_integer,0,mpi_any_tag,mpi_comm_world,&
             &status,ierr)
     end if

     print *, ntasks, id, i

  end do

  call mpi_finalize(ierr)

end program probably_works
