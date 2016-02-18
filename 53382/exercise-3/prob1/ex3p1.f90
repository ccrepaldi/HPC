program ackermann
  implicit none
  integer :: m,n,ack

  write(6,*) "Insert a value for the integer m:"
  read(5,*) m
  write(6,*) "Insert a value for the integer n:"
  read(5,*) n

  write(*,*) "A(m,n) =", ack(m,n)

end program ackermann

recursive function ack(m,n) result(a)
  integer, intent(in) :: m,n ! variable can enter but cannot be changed
  integer :: a
  if (m==0) then
     a=n+1
  else if (n==0) then
     a=ack(m-1,1)
  else
     a=ack(m-1,ack(m,n-1))
  end if
end function ack
