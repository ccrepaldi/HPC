
!---------------------------------!
! Compilation by Intel Fortran:   !
! ifort -openmp omp_hello.f90     !
!---------------------------------!
! Compilation by GNU Fortran:     !
! gfortran -fopenmp omp_hello.f90 !
!---------------------------------!

program omp_hello
  !$ use omp_lib
  implicit none
  integer :: tid,threads,i,n

  tid=-1
  threads=-1
  n=0

  !$ print *,'Procs = ',omp_get_num_procs()
  !!!!!$ call omp_set_num_threads(4)

  !$omp parallel private(tid),shared(i),shared(n)
  !$ tid = omp_get_thread_num()
  !$ threads=omp_get_num_threads()
  n=n+1
  i=tid
  print *, 'Thread = ', tid,' out of ',threads,i,n
  !$omp end parallel  

  print *,'I''m the master:',tid,i

end program omp_hello
