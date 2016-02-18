program optimization1
  implicit none
  integer :: i
  integer, parameter :: n = 10000000
  real :: a(n), b(n),t1,t2
                                      
  a(:)=42.0
  b(:)=42.0
  call cpu_time(t1)
  do i=1,n-1
     if (i<500) then
        a(i)=4.0*b(i)+b(i+1)
     else
        a(i)=4.0*b(i+1)+b(i)
     end if
  end do
  call cpu_time(t2)
  print *, "Elapsed time:", t2-t1
end program optimization1
