program unroll_test
  implicit none
  integer :: m,n=10000,i
  integer(kind=4), dimension(10000) :: a
  real(kind=16) :: t1,t2

  ! Unroll: depth=1
  call cpu_time(t1)
  m = mod(n,1)
  do i=1,m
     a(i) = i
  end do

  do i = m+1,n,1
     a(i+0) = i+0
  end do
  call cpu_time(t2)
  print *, "Elapsed time:", t2-t1
  ! Unroll: depth=2
  call cpu_time(t1)
  m = mod(n,2)
  do i=1,m
     a(i) = i
  end do

  do i = m+1,n,2
     a(i+0) = i+0
     a(i+1) = i+1
  end do
  call cpu_time(t2)
  print *, "Elapsed time:", t2-t1
  ! Unroll: depth=4
  call cpu_time(t1)
  m = mod(n,4)
  do i=1,m
     a(i) = i
  end do

  do i = m+1,n,4
     a(i+0) = i+0
     a(i+1) = i+1
     a(i+2) = i+2
     a(i+3) = i+3
  end do
  call cpu_time(t2)
  print *, "Elapsed time:", t2-t1
  ! Unroll: depth=8
  call cpu_time(t1)
  m = mod(n,8)
  do i=1,m
     a(i) = i
  end do

  do i = m+1,n,8
     a(i+0) = i+0
     a(i+1) = i+1
     a(i+2) = i+2
     a(i+3) = i+3
     a(i+4) = i+4
     a(i+5) = i+5
     a(i+6) = i+6
     a(i+7) = i+7
  end do
  call cpu_time(t2)
  print *, "Elapsed time:", t2-t1
  ! Unroll: depth=16
  call cpu_time(t1)
  m = mod(n,16)
  do i=1,m
     a(i) = i
  end do

  do i = m+1,n,16
     a(i+0) = i+0
     a(i+1) = i+1
     a(i+2) = i+2
     a(i+3) = i+3
     a(i+4) = i+4
     a(i+5) = i+5
     a(i+6) = i+6
     a(i+7) = i+7
     a(i+8) = i+8
     a(i+9) = i+9
     a(i+10) = i+10
     a(i+11) = i+11
     a(i+12) = i+12
     a(i+13) = i+13
     a(i+14) = i+14
     a(i+15) = i+15
  end do
  call cpu_time(t2)
  print *, "Elapsed time:", t2-t1
  ! Unroll: depth=32
  call cpu_time(t1)
  m = mod(n,32)
  do i=1,m
     a(i) = i
  end do

  do i = m+1,n,32
     a(i+0) = i+0
     a(i+1) = i+1
     a(i+2) = i+2
     a(i+3) = i+3
     a(i+4) = i+4
     a(i+5) = i+5
     a(i+6) = i+6
     a(i+7) = i+7
     a(i+8) = i+8
     a(i+9) = i+9
     a(i+10) = i+10
     a(i+11) = i+11
     a(i+12) = i+12
     a(i+13) = i+13
     a(i+14) = i+14
     a(i+15) = i+15
     a(i+16) = i+16
     a(i+17) = i+17
     a(i+18) = i+18
     a(i+19) = i+19
     a(i+20) = i+20
     a(i+21) = i+21
     a(i+22) = i+22
     a(i+23) = i+23
     a(i+24) = i+24
     a(i+25) = i+25
     a(i+26) = i+26
     a(i+27) = i+27
     a(i+28) = i+28
     a(i+29) = i+29
     a(i+30) = i+30
     a(i+31) = i+31
  end do
  call cpu_time(t2)
  print *, "Elapsed time:", t2-t1
  
end program unroll_test
