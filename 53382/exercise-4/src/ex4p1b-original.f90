program optimization2
  implicit none
  integer, parameter :: n=10000
  real :: a(n,n),b(n,n),c(n),t1,t2,tt1,tt2
  integer :: i,j,m
  b(:,:)=24.0
  c(:)=2.0

  call cpu_time(t1)
  do i=1,n
     do j=1,n
        a(i,j)=b(i,j)/c(i)
     end do
  end do
  call cpu_time(t2)

  print *, "Elapsed time:", t2-t1

end program optimization2
