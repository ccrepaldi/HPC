program optimization1
  implicit none
  integer :: i,m
  integer, parameter :: n = 10000000
  real :: a(n), b(n),t1,t2
  
  a(:)=42.0
  b(:)=42.0
  call cpu_time(t1)
  
  ! Unroll: depth=8
  m = mod(499,8)
  do i=1,m
     a(i) = 4.0*b(i)+b(i+1)
  end do

  do i = m+1,499,8
     a(i) = 4.0*b(i+0)+b(i+0+1)
     a(i+1) = 4.0*b(i+1)+b(i+1+1)
     a(i+2) = 4.0*b(i+2)+b(i+2+1)
     a(i+3) = 4.0*b(i+3)+b(i+3+1)
     a(i+4) = 4.0*b(i+4)+b(i+4+1)
     a(i+5) = 4.0*b(i+5)+b(i+5+1)
     a(i+6) = 4.0*b(i+6)+b(i+6+1)
     a(i+7) = 4.0*b(i+7)+b(i+7+1)
  end do

  ! Unroll: depth=8
  m = mod(n-1,8)
  do i=500,m
     a(i) = 4.0*b(i+1)+b(i)
  end do

  do i = m+1,n-1,8
     a(i+0)=4.0*b(i+0+1)+b(i+0)
     a(i+1)=4.0*b(i+1+1)+b(i+1)
     a(i+2)=4.0*b(i+2+1)+b(i+2)
     a(i+3)=4.0*b(i+3+1)+b(i+3)
     a(i+4)=4.0*b(i+4+1)+b(i+4)
     a(i+5)=4.0*b(i+5+1)+b(i+5)
     a(i+6)=4.0*b(i+6+1)+b(i+6)
     a(i+7)=4.0*b(i+7+1)+b(i+7)
  end do

!  a(1:499) = 4.0*b(1:499)+b(2:500)
!  a(500:n) = 4.0*b(501:n+1)+b(500:n)
  call cpu_time(t2)
  print *, "Elapsed time:", t2-t1
end program optimization1
