program simpleprofiling
  integer,parameter :: rk=selected_real_kind(5,100)
  integer,parameter :: NMAX=100000,MMAX=1000
  real(rk),parameter :: dx=0.00001
  real(rk) :: sum1,sum2
  integer :: i,j

  sum1=0.0
  sum2=0.0

  do j=0,MMAX
     do i=0,NMAX
        sum1=sum1+dx*i*j
     end do
  end do
     
  do j=0,MMAX
     do i=0,NMAX
        sum2=sum2+cos(exp(sin(dx*i*j)))
     end do
  end do

  print *,sum1,sum2

end program simpleprofiling
