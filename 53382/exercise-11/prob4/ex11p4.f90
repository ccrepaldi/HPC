program omp_pi
  !$ use omp_lib
  implicit none
  integer,parameter :: rk=selected_real_kind(10,40)    ! double precision
  integer,parameter :: n=100, m=100, p=100
  integer :: tid,threads,i,j,k,chunk=10
  real(rk) :: a(n,m), b(m,p), c(n,p),sum,wt0,wt1

  tid=-1
  threads=-1
  
  ! initializing matrices

  a(:,:)=42.0
  b(:,:)=66.7
  
  !$ call omp_set_num_threads(10)
  !$omp parallel private(tid,i,j,k),shared(a,b,c,threads,chunk)
  
  !$ tid = omp_get_thread_num()
  !$ threads=omp_get_num_threads()

  !$ write (6,'(a,2x,i0,2x,a,i0)') trim("Number of threads:"), threads, trim(" ---> I am thread #"), tid

  !$ wt0=omp_get_wtime()
  !$omp parallel do schedule(static,chunk)
  do i=1,n
     !$ write (6,'(a, 1x, i0, a, i0)') trim("Thread"), tid, trim(" is taking care of row #"), i
     do j=1,p
        sum=0.0
        do k=1,m
           sum=sum+(a(i,k)*b(k,j))
        end do
        c(i,j)=sum
     end do
  end do
  !$omp end parallel
  !$ wt1=omp_get_wtime()
  !$ wt1=wt1-wt0

  write (6,'(a,2x,f10.5,a)') trim("Wall clock time = "), wt1, trim("s")

 ! do i=1,n
 !    write (6,*) [(c(i,j), j=1,p)]
 ! end do

end program omp_pi

