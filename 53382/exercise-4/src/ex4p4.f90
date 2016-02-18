program vectorize_test
  implicit none
  integer, parameter :: n=10000
  real :: a(n,n),b(n,n),c(n),d=255,t1,t2
  integer :: i,j

  call cpu_time(t1)
  do i=1,n
     b(i,:)=(i*20)/17.0
     c(i)=i*23
  end do
  

  do i=1,n
     do j=1,n
        a(i,j)=c(i)*b(j,i)+d
     end do
  end do
  call cpu_time(t2)

  print *, "Elapsed time:", t2-t1

end program vectorize_test
