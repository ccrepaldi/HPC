program lineardouble

! solving the matrix equation A*x=b using LAPACK
  implicit none

  integer,parameter :: rk=8
  real(rk),allocatable :: a(:,:),b(:),x(:),a0(:,:),b0(:),b1(:)
  integer :: n,i,j,ok
  integer,allocatable :: pivot(:)

  read(5,*) n
  allocate(a(n,n),b(n),x(n),pivot(n))
  allocate(a0(n,n),b0(n),b1(n))
  do i=1,n
     read(5,*) (a(i,j),j=1,n)
  enddo
  read(5,*) (b(i),i=1,n)
  a0=a
  b0=b


  write(6,'(/a)') 'Before'
  write(6,'(/a)') 'A'
  do i=1,n
     write(6,'(20g16.8)') (a(i,j),j=1,n)
  enddo
  write(6,'(/a)') 'b'
  write(6,'(20g16.8)') b

  call dgesv(n, 1, a, n, pivot, b, n, ok)
  x=b

  write(6,'(/a)') 'After'
  write(6,'(/a)') 'A'
  do i=1,n
     write(6,'(20g16.8)') (a(i,j),j=1,n)
  enddo
  write(6,'(/a)') 'x'
  write(6,'(20g16.8)') x

  write(6,'(//)')
  write(6,'(a,20g16.8)') 'A*x      ',matmul(a0,x)
  write(6,'(a,20g16.8)') 'b        ',b0
  write(6,'(a,20g16.6)') 'Residual ',b0-matmul(a0,x)


end program lineardouble

