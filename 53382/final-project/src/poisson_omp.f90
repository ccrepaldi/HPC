module sizes
  integer,parameter :: rk0=selected_real_kind(5,10)    ! single precision
  integer,parameter :: rk=selected_real_kind(10,40)    ! double precision
  integer,parameter :: rk2=selected_real_kind(30,200)  ! quadruple precision
  integer,parameter :: ik=selected_int_kind(12)        ! long int
  integer,parameter :: MAXBUF=200
end module sizes

module gfunctionmod
  use sizes
  implicit none

contains

  real(rk) function g_function(x,y) ! this is the function g(x,y)
    implicit none
    real(rk), intent(in) :: x,y

    g_function=0.0

  end function g_function

end module gfunctionmod

program poisson_parallel
  !$ use omp_lib
  use sizes
  use gfunctionmod
  implicit none
  integer,parameter :: N=128 ! system size
  real(rk) :: gamma=1.9378666 ! over relaxation parameter
  real(rk),dimension(N,N) :: f,g,f_old ! f (old and new) and g matrices
  integer :: i,j ! indeces
  integer,parameter :: maxiter=50000 ! maximum number of iterations
  integer :: iter ! current iteraction
  real(rk) :: temp ! temporary variable with the difference between the old and
  ! and new value of an element of the matrix
  !
  real(rk) :: diff ! biggest difference between the old and new values of the
  ! matrix elements
  integer :: narg ! number of command line arguments
  integer :: threads,nt ! openmp variables
  character(len=80) :: arg ! command line argument
  real(rk) :: wt0,wt1 ! wallclock times

  narg=command_argument_count()

  !$ if (narg==1) then
  !$    call get_command_argument(1,arg); read(arg,*) nt
  !$    if (nt>0) call omp_set_num_threads(nt) 
  !$ else
  !$ call omp_set_num_threads(4)
  !$ end if
  !$ print '(a,x,i0)', "CPUs (cores) available:", omp_get_num_procs()
  !$omp parallel
  !$omp master
  !$ print '(a,x,i0,x,a)', "Running with", omp_get_num_threads(), "threads."
  !$omp end master
  !$omp end parallel
  !$ print "(50('-'))"

  print '(a,i0)', "N = ", N
  print '(a,f7.3)', "gamma = ", gamma

  print "(50('-'))"

  !$ wt0=omp_get_wtime()

  ! lets fill the g matrix first

  do i=1,N
     !$omp parallel do private(j)
     do j=1,N
        g(i,j)=g_function((1.0_rk*i)/N,(1.0_rk*j)/N)
     end do
     !$omp end parallel do
  end do

  ! lets fill f with null values in order to remove the
  ! trash from the array address

  do i=1,N
     !$omp parallel do private(j)
     do j=1,N
        f(i,j)=0.0
     end do
     !$omp end parallel do
  end do

  ! let's update f with some border conditions

  !$omp parallel do private(i)
  do i=1,N
     f(i,1)=-(i-(N/2.0))**2/(N/2.0)**2 +1.0 ! parabola
  end do
  !$omp end parallel do

  !$omp parallel do private(i)
  do i=1,N
     f(i,N)=0.0
  end do
  !$omp end parallel do

  !$omp parallel do private(j)
  do j=1,N
     f(1,j)=0.0
  end do
  !$omp end parallel do

  !$omp parallel do private(j)
  do j=1,N
     f(N,j)=0.0
  end do
  !$omp end parallel do

  ! Red-black algorithm ************************************

  ! First we update the red sublattice values for the first time
  ! Note that we do this update without touching the values at the borders
  ! because of the border conditions

  !$omp parallel do private(i,j)
  do i=2,N-1
     do j=2,N-2
        if (mod(i+j,2)/=0) then ! we have a red dot
           f(i,j)= 0.0 ! some guess value
        end if
     end do
  end do
  !$omp end parallel do

  iter=0

  coverge_loop: do

     diff=0.0 ! biggest difference between old and new values

     !$omp parallel do private(i,j)
     do i=1,N
        do j=1,N
           f_old(i,j)=f(i,j) ! save the old version of f
        end do
     end do
     !$omp end parallel do

     !$omp parallel do private(i,j,temp) shared(diff,f)
     do i=2,N-1
        do j=2,N-1
           if (mod(i+j,2)==0) then ! we have a black dot
              f(i,j)=(1.0_rk-gamma)*f_old(i,j)+(gamma/4.0_rk)*&
                   &(f(i+1,j)+f(i-1,j)+f(i,j+1)+f(i,j-1))&
                   &-(gamma/(4.0_rk*N**2)*g(i,j))
              temp=abs(f_old(i,j)-f(i,j))
              if (temp > diff) then
                 diff=temp
              end if
           end if
        end do
     end do
     !$omp end parallel do

     !$omp parallel do private(i,j,temp) shared(diff,f)
     do i=2,N-1
        do j=2,N-1
           if (mod(i+j,2)/=0) then ! we have a red dot
              f(i,j)=(1.0_rk-gamma)*f_old(i,j)+(gamma/4.0_rk)*&
                   &(f(i+1,j)+f(i-1,j)+f(i,j+1)+f(i,j-1))&
                   &-(gamma/(4.0_rk*N**2)*g(i,j))
              temp=abs(f_old(i,j)-f(i,j))
              if (temp > diff) then
                 diff=temp
              end if
           end if
        end do
     end do
     !$omp end parallel do

     iter=iter+1

     if (diff<=0.000001) then ! the method converged
        print '(a)', trim('Method converged!')
        print '(a,x, i0)', trim('Iter ='), iter
        print '(a,x,f10.6)', trim('Diff ='), diff
        exit
     end if

     if (iter==maxiter) then
        print *, "The method reached its maximum iteration number. Iter = ", iter
        exit
     end if

  end do coverge_loop

  !$ wt1=omp_get_wtime()
  !$ print '(a,f10.5)', "Wall clock time: ",wt1-wt0

end program poisson_parallel
