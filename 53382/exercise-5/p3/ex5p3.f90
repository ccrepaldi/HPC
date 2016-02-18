program command_line
  implicit none
  integer :: i
  real :: t1,t2
  integer, parameter :: sp = selected_real_kind(6, 37), n=10000000
  integer :: count1,count2, count_rate, count_max
  real(kind=sp), dimension(n) :: a

  a(:)=42.0

  open(2,file='un1.txt',form='unformatted',access='stream',status='replace')
  call system_clock(count1, count_rate)
  do i=1,n
     write(2) a(i)
  end do
  call system_clock(count2, count_rate)
  write(6,*) count2-count1
  close(unit=1,status='keep')

  open(2,file='un2.txt',form='formatted',access='stream',status='replace')
  call system_clock(count1, count_rate)
  do i=1,n
     write(2,'(g0.10)') a(i)
  end do
  call system_clock(count2, count_rate)
  write(6,*) count2-count1
  close(unit=1,status='keep')
  
end program command_line
