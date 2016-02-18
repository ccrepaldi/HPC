program ex3p3
  implicit none
  integer,parameter :: N=15000
  integer :: a(N,N)
  integer :: i,j
  real(kind=10) :: t1,t2

  call cpu_time(t1)! Begin measurement
  do j=1,N
     do i=1,N
        a(i,j)=(i+j)/2
     end do
  end do
  call cpu_time(t2)! End measurement
  print *, "Time spent =", t2-t1
  !print '(10i8)',a(1:N:10000,1:N:10000)

end program ex3p3

