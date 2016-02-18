
program eigendouble

! Eigenvalues of a square matrix using LAPACK DGEEV

  implicit none

  
  integer,parameter :: rk=selected_real_kind(10,40)
  integer :: n,i,j,ok,lda,ldvl,ldvr,lwork
  real(rk),allocatable :: a(:,:),wr(:),wi(:),vl(:,:),vr(:,:),work(:)
  logical :: symmetric,uppertriangular

  symmetric=.false.
  uppertriangular=.false.

  read(5,*) n
  lwork=4*n
  lda=n
  ldvl=n
  ldvr=n
  !allocate(a(n,n),wr(n),wi(n),vl(n,n),vr(n,n),work(lwork))
  allocate(a(n,n),wr(n),wi(n),vl(n,n),vr(1,1),work(lwork))
  do i=1,n
     read(5,*) (a(i,j),j=1,n)
  enddo

  if (symmetric) a=(a+transpose(a))/2.0
  if (uppertriangular) forall (i=1:n,j=1:n,i>j) a(i,j)=0.0

  write(6,'(/a)') 'A'
  do i=1,n
     write(6,'(20g16.8)') (a(i,j),j=1,n)
  enddo

  call dgeev('V','N',n,a,n,wr,wi,vl,n,vr,n,work,lwork,ok)
  write(6,'(/,a)') 'Eigenvalues of A:'
  write(6,'(2g16.8)') (wr(i),wi(i),i=1,n)
  write(6,'(/a)') 'Eigenvectors of A'
  do i=1,n
     write(6,'(20g16.8)') (vl(i,j),j=1,n)
  enddo


end program eigendouble

