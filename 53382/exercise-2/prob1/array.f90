program array
  implicit none
  integer :: i, it, a(10)

  write(6,*) "Please give the number iteractions:"
  read(*,*) it

  do i=1,it
     print *, i, a(i)
  end do

end program array
