program command_line
  implicit none
  integer :: i
  real :: t1,t2
  integer, parameter :: sp = selected_real_kind(6, 37), n=10000000
  real(kind=sp), dimension(n) :: a

  a(:)=42.0

  open(2,file='un1.txt',form='unformatted',access='stream',status='replace')

  call cpu_time(t1)
  do i=1,n
     write(2) a(i)
  end do
  call cpu_time(t2)

  write(6,*) "Elapsed time:", t2-t1

  close(unit=1,status='keep')

  open(2,file='un2.txt',form='formatted',access='stream',status='replace')
  
  call cpu_time(t1)
  do i=1,n
     write(2,'(g0.10)') a(i)
  end do
  call cpu_time(t2)

  write(6,*) "Elapsed time:", t2-t1
  
end program command_line
