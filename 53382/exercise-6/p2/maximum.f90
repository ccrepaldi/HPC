program maximum
  implicit none
  integer,parameter :: NMAX=20000000,M=20,K=2001,rk=8
  real(rk) :: a(NMAX),x,t1,t2
  integer :: i

  forall (i=1:NMAX)
     a(i)=exp(real(mod((i-NMAX),M),rk)/(10.0*M))+sin(real(mod(i,1000*M)/K,rk))
  end forall

  x=-huge(x)
  call cpu_time(t1)

  do i=1,NMAX
     if (a(i)>=x) x=a(i)
  end do

  call cpu_time(t2)

  print *,t2-t1,x

end program maximum
