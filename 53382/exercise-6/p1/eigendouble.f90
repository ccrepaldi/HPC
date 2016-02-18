
program eigendouble

! Eigenvalues of a square matrix using LAPACK DGEEV

  implicit none

  
  integer,parameter :: rk=selected_real_kind(10,40)
  integer :: n,i,j,ok,lda,ldvl,ldvr,lwork
  real(rk),allocatable :: a(:,:),wr(:),wi(:),vl(:,:),vr(:,:),work(:)
  real :: tf,ti
  character(len=80) :: arg
  integer :: seed

  ! --- Command line arguments: matrix size, RNG seed

  if (command_argument_count()/=2) then
     call get_command_argument(0,arg)
     write(0,'(a)') 'usage: '//trim(arg)//' n seed'
     stop
  end if
  call get_command_argument(1,arg); read(arg,*) n
  call get_command_argument(2,arg); read(arg,*) seed


  ! --- Allocate space

  lwork=4*n
  lda=n
  ldvl=n
  ldvr=n
  allocate(a(n,n),wr(n),wi(n),vl(n,n),vr(n,n),work(lwork))

  ! --- Fill the matrix with randomn numbers

  call fill_matrix()

  ! ---- The eigenvalue calculation proper ----

  call cpu_time(ti)
  call dgeev('V','V',n,a,n,wr,wi,vl,n,vr,n,work,lwork,ok)
  call cpu_time(tf)

  write(6,*) tf-ti, n
  !write(6,'(50("-"))')

  ! --- Print eigenvalues to file "evalues.datf"

  open(100,file="evalues.datf",status="replace")
  write(100,'(2g16.8)') (wr(i),wi(i),i=1,n)
  close(100)
  

contains

  subroutine fill_matrix()
    implicit none
    integer :: i,j
    real :: x
    x=ran(seed) ! This is only for setting the seed
    do i=1,n
       do j=1,n
          a(i,j)=ran(0)
       end do
    end do
    return
  end subroutine fill_matrix

  ! --- A simple random number generator

  real(rk) function ran(seed)
    implicit none
    integer, intent(in) :: seed
    integer,save :: ix
    integer,parameter :: ii =127773, jj=2147483647
    integer,parameter :: i1=16807, i2=2836
    real(rk),parameter :: p = 4.656612875d-10
    integer :: k1
    if (seed/=-1.and.seed/=0) ix=seed
    k1 = ix/ii
    ix = i1*(ix-k1*ii)-k1*i2
    if (ix<0) ix = ix+jj
    ran = ix*p
    return
  end function ran
  

end program eigendouble

