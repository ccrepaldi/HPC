program matrixproducts
  implicit none
  real :: a(1000,1000),b(1000,1000),c(1000,1000)
  integer :: i,j

  a(:,:)=42.0
  b(:,:)=a(:,:)

  
  call product1(a,b,c)
  
  call product2(a,b,c)
  
  call product3(a,b,c)
 
end program matrixproducts

subroutine product1(a,b,c)
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)
  real :: a(1000,1000),b(1000,1000),c(1000,1000)
  real(kind=dp) :: ti,tf
  integer :: i,j,k
  call cpu_time(ti)
  c(:,:)=0.0
  do i=1,1000
     do j=1,1000
        do k=1,1000
           c(i,j)=c(i,j)+a(i,k)*b(k,j)
        end do
     end do
  end do
  call cpu_time(tf)
  print *, "time:", tf-ti
end subroutine product1

subroutine product2(a,b,c)
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)
  real :: a(1000,1000),b(1000,1000),c(1000,1000)
  real(kind=dp) :: ti,tf
  integer :: i,j,k
  call cpu_time(ti)
  c(:,:)=0.0
  do i=1,1000
     do j=1,1000
        do k=1,1000
           c(i,j)=c(i,j)+a(j,k)*b(k,j)
        end do
     end do
  end do
  call cpu_time(tf)
  print *, "time:", tf-ti
end subroutine product2

subroutine product3(a,b,c)
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)
  real :: a(1000,1000),b(1000,1000),c(1000,1000)
  real(kind=dp) :: ti,tf
  integer :: i,j,k
  call cpu_time(ti)
  c(:,:)=0.0
  do i=1,1000
     do j=1,1000
        do k=1,1000
           c(i,j)=c(i,j)+a(i,k)*b(k,i)
        end do
     end do
  end do
  call cpu_time(tf)
  print *, "time:", tf-ti
end subroutine product3
