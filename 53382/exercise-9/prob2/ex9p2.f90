! The program only works if, in this case, the number of processes is 8
! Otherwise it will just mess with the matrix elements
! See more info in the README file
program mpi_transpose
  use mpi
  implicit none
  integer, parameter :: ROOT=0
  character(len=80) :: a
  integer :: id,ntasks,rc,n,i,j,ncol
  integer, allocatable :: matrix(:,:), piece0(:), piece1(:), lmatrix(:)

  call mpi_init(rc)
  call mpi_comm_size(mpi_comm_world,ntasks,rc)
  call mpi_comm_rank(mpi_comm_world,id,rc)

  n=8 ! matrix shall be 8x8
  ncol=8/ntasks ! number of columns that shall be stored in each process
  allocate(matrix(n,n),piece0(n*ncol),piece1(n*ncol),lmatrix(n*n))

  ! first we create the matrix
  do i=1,n
     do j=1,n
        matrix(i,j)=10*i+j
     end do
  end do

  ! now we print the original matrix to check on later
  if (id==ROOT) then
     do i=1,n
        print '(8(i3))', matrix(i,:) ! Matrix OK
     end do
     print '(A)', new_line(a)

     lmatrix=reshape(matrix,(/ n*n /))

    ! print '(i0)', lmatrix ! Linear matrix OK
     
  end if

  call mpi_scatter(lmatrix,n*ncol,mpi_integer,piece0,n*ncol,mpi_integer,ROOT,mpi_comm_world,rc)

 ! if (id==ntasks-2) print *, id, new_line(a), piece0 ! Scattering OK

  call mpi_alltoall(piece0,(n*ncol)/ntasks,mpi_integer,piece1,(n*ncol)/ntasks,mpi_integer,mpi_comm_world,rc)

 ! if (id==ntasks-2) print *, id, new_line(a), piece1 ! All2All not ok -> it is not tranposing the matrix columns

  call mpi_gather(piece1,n*ncol,mpi_integer,lmatrix,n*ncol,mpi_integer,ntasks-1,mpi_comm_world,rc) ! Gather OK

  if (id==ntasks-1) then
    ! print '(i0)', lmatrix ! OK 
     matrix=reshape(lmatrix,(/ n, n /),order=(/2,1/))
     print '(8(i3))', matrix ! Only OK if number of tasks = ncol or nrows of the matrix
  end if

  call mpi_finalize(rc)

end program mpi_transpose
