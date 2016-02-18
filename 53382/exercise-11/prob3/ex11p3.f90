program omp_pi
  use omp_lib
  implicit none
  integer,parameter :: rk=selected_real_kind(10,40)    ! double precision
  integer :: tid,threads,i,n,nt
  real(rk) :: pi
  character(len=80) :: argu
  ! get the command line arguments
  call getarg(1,argu); read(argu,*) nt ! number of threads
  call getarg(2,argu); read(argu,*) n ! number of terms
  ! set the desired number of threads
  call omp_set_num_threads(nt)
  ! initialize the pi value before sum as 0
  pi=0.0

  ! begin parallel sum

  !$omp parallel do reduction(+:pi)
  do i=0,n
     pi = pi + (-1.0)**i/(2.0*i+1)
  end do
  !$omp end parallel do

  ! end parallel sum
  
  pi=pi*4

  ! prints the number of threads and the pi value
  print *, nt, pi

end program omp_pi

