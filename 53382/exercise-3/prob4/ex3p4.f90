program precision_speed
  implicit none
  integer :: n

  write(*,*) "Please insert a value for the size of the array:"
  read(*,*) n

  call arrays(n)

end program precision_speed

subroutine arrays (n)
  implicit none
  integer :: i,j
  integer, intent(in) :: n
  integer, parameter :: sp = selected_real_kind(5,10) ! single precision
  integer, parameter :: dp = selected_real_kind(10,40) ! double precision
  integer, parameter :: qp = selected_real_kind(20,1000) ! quadruple precision
  real(sp), dimension (n) :: a
  real(dp), dimension(n) :: b
  real(qp), dimension(n) :: c
  real :: sum
  real :: t1,t2,t3,t4,t5,t6,tcpu

  call cpu_time(t1)

  sum = 0.0
  do i=1,(n-1)
     do j=(i+1),n
        sum = sum + (a(j)-a(i))
     end do
  end do

  call cpu_time(t2)
  tcpu = t2-t1
  write(*,*) "CPUtime for single-precision array =", tcpu

  call cpu_time(t3)
  sum =0.0
  do i=1,(n-1)
     do j=(i+1),n
        sum = sum + (b(j)-b(i))
     end do
       end do
  call cpu_time(t4)
  tcpu = t4-t3
  write(*,*) "CPUtime for double-precision array =", tcpu

  call cpu_time(t5)
  sum =0.0
  do i=1,(n-1)
     do j=(1+1),n
        sum = sum + (c(j)-c(i))
     end do
       end do
  call cpu_time(t6)
  tcpu = t6-t5
  write(*,*) "CPUtime for quadruple-precision array =", tcpu

end subroutine arrays
