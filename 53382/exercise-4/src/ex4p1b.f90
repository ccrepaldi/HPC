program optimization2
  implicit none
  integer, parameter :: n=10000
  real :: a(n,n),b(n,n),c(n),d,t1,t2
  integer :: i,j
  b(:,:)=23.0
  c(:)=7.0

  call cpu_time(t1)
  a(:,:)=b(:,:)
  do i=1,n,8
     a(i,:)=a(i,:)/c(i)
     a(i+1,:)=a(i+1,:)/c(i+1)
     a(i+2,:)=a(i+2,:)/c(i+2)
     a(i+3,:)=a(i+3,:)/c(i+3)
     a(i+4,:)=a(i+4,:)/c(i+4)
     a(i+5,:)=a(i+5,:)/c(i+5)
     a(i+6,:)=a(i+6,:)/c(i+6)
     a(i+7,:)=a(i+7,:)/c(i+7)
  end do 
  call cpu_time(t2)

  print *, "Elapsed time:", t2-t1

end program optimization2
