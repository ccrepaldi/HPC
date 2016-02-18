program omp_pardo
  use omp_lib
  implicit none
  integer,parameter :: rk=selected_real_kind(10,40)    ! double precision
  integer :: tid,threads,i,n,nt
  real(rk) :: s,wt0,wt1
  character(len=80) :: argu

  call getarg(1,argu); read(argu,*) nt
  call getarg(2,argu); read(argu,*) n

  call omp_set_num_threads(nt)

  s=0.0
  wt0=omp_get_wtime()  
  !$omp parallel do reduction(+:s)
  do i=1,n
     s=s+log(real(i,rk))
  end do
  !$omp end parallel do
  wt1=omp_get_wtime()
  wt1=wt1-wt0
  print '(i0,2x,2g20.10)',nt,wt1,s

end program omp_pardo

