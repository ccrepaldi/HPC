! ************************************************************************
! Program MD1D with MPI implementation
! Tools for HPC (Autumn 2015) - Helsingin Yliopisto
! Student: Caike Crepaldi <caike.crepaldi@helsinki.fi>
! Program compiled and ran in pangolin server, root directory
! ************************************************************************
module sizes
  integer,parameter :: rk0=selected_real_kind(5,10)    ! single precision
  integer,parameter :: rk=selected_real_kind(10,40)    ! double precision
  integer,parameter :: rk2=selected_real_kind(30,200)  ! quadruple precision
  integer,parameter :: ik=selected_int_kind(12)        ! long int
  integer,parameter :: MAXBUF=200
end module sizes

module mpi_var ! module with mpi only variables and routines
  use mpi
  use sizes
  integer,parameter :: ROOT=0,tag=50
  integer :: myid,ntasks,source_id,dest_id,ierr
  integer,dimension(mpi_status_size) :: status
  real(rk),dimension(2) :: msg

contains
  ! subroutine that allows the neighboring processes to comunicate with each other
  subroutine neighboor_communication (lim1,lim2,myid,ntasks,msg)
    implicit none
    integer, intent(in) :: myid, ntasks
    real(rk), intent(in) :: lim1, lim2 ! lim1 = x(1) and lim2 = x(niter)
    real(rk), intent(out) :: msg(2) ! msg(1) = x(0) and msg(2) = x(niter+1)
    
    ! source ------> dest
    dest_id=myid+1
    if (dest_id.gt.ntasks-1) dest_id=ROOT ! (ntasks-1) ------> ROOT
    source_id=myid-1
    if (source_id.lt.0) source_id=(ntasks-1) ! ROOT <------ (ntasks-1)
    
    if (mod(myid,2)==0) then ! trying to avoid deadlocks
       call mpi_send(lim2,1,mpi_double_precision,dest_id,tag,mpi_comm_world,ierr)
       call mpi_recv(msg(1),1,mpi_double_precision,source_id,&
            &tag,mpi_comm_world,status,ierr)
       call mpi_send(lim1,1,mpi_double_precision,source_id,tag,mpi_comm_world,ierr)
       call mpi_recv(msg(2),1,mpi_double_precision,dest_id,&
            &tag,mpi_comm_world,status,ierr)
    else
       call mpi_recv(msg(1),1,mpi_double_precision,source_id,&
            &tag,mpi_comm_world,status,ierr)
       call mpi_send(lim2,1,mpi_double_precision,dest_id,tag,mpi_comm_world,ierr)
       call mpi_recv(msg(2),1,mpi_double_precision,dest_id,&
            &tag,mpi_comm_world,status,ierr)
       call mpi_send(lim1,1,mpi_double_precision,source_id,tag,mpi_comm_world,ierr)
    end if

  end subroutine neighboor_communication

end module mpi_var

program md1d
  use mpi
  use mpi_var
  use sizes
  implicit none  

  real(rk),allocatable :: x(:)   !    atom positions
  real(rk),allocatable :: v(:)   !    velocities
  real(rk),allocatable :: v0(:)  !    previous veloocities (leap frog needs them)
  real(rk),allocatable :: a(:)   !    accelerations
  real(rk),allocatable :: ep(:)  !    potential energies
  real(rk),allocatable :: ek(:)  !    kinetic energies
  real(rk),allocatable :: allcoord(:), allep(:) ! all coordinates and poten. energies

  real(rk) :: epsum,eksum        ! system energies

  real(rk) :: dt                 ! time step 
  real(rk) :: vsc                ! mean initial velocity
  real(rk) :: box                ! system size
  integer :: nat                 ! number of atoms
  integer :: maxt                ! number of time steps simulated
  integer :: eout                ! energy output interval
  integer :: cout                ! coordinate output interval (lot of data, beware!)
  
  integer :: i,n,ia
  character(len=MAXBUF) :: arg
  integer(ik) :: seed=987234

  integer :: niter               ! number of iterations per processor
  integer :: begindex            ! index of the first element in the processor
  integer :: endindex            ! index of the last element in the processor
  real(rk) :: vsum               ! sum of the velocities of all atoms
  real(rk) :: t0, t1, tused      ! wall clock times
  logical :: time_print ! if it is .true. the program prints the execution time at the end

  ! Get number of atoms, time step and simulation length from command line
  ia=command_argument_count()
  if (ia<4.or.ia>6) then
     if (myid==ROOT) then
        call get_command_argument(0,arg)
        write(6,'(/,a,a,a)') 'usage: ',trim(arg),' nat dt maxt vsc [eout [cout]]'
        write(6,'(a)')       '    nat  = number of atoms'
        write(6,'(a)')       '    dt   = time step'
        write(6,'(a)')       '    maxt = number of time steps in simulation'
        write(6,'(a)')       '    vsc  = mean velocity of atoms in the beginning (''temperature'')'
        write(6,'(a)')       '    eout = interval for printing energies to stdout'
        write(6,'(a,/)')     '    cout = interval for printing coordinates to ''fort.10'''
     end if
     stop
  end if
  cout=0
  eout=0 ! eout=1 ! changed in order to print only the times
  call get_command_argument(1,arg); read(arg,*) nat
  call get_command_argument(2,arg); read(arg,*) dt
  call get_command_argument(3,arg); read(arg,*) maxt
  call get_command_argument(4,arg); read(arg,*) vsc
  if (ia==4) then ! if the cout and eout are not passed as arguments
     time_print=.true. ! print the execution time
  else
     time_print=.false.
  end if
  if (ia>4) then
     call get_command_argument(5,arg); read(arg,*) eout
  end if
  if (ia>5) then
     call get_command_argument(6,arg); read(arg,*) cout
  end if

  ! begin mpi initialization
  call mpi_init(ierr)
  if (ierr/=mpi_success) then
     print *,'MPI initialization failed.'
     stop
  end if
  call mpi_comm_size(mpi_comm_world,ntasks,ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  ! end mpi initialization

  if (mod(nat,ntasks)/=0) then
     if (myid==ROOT) print *, 'Number of atoms must be exactly divisible by the number of processes.'
     call mpi_finalize(ierr)
     stop
  end if

  t0=MPI_Wtime()

  niter = nat/ntasks
  begindex = myid*niter+1
  endindex = niter+begindex-1
  
  ! allocate arrays in each process with size equals to nat/ntasks
  allocate(x(niter),v(niter),v0(niter),a(niter),ep(niter),ek(niter),allcoord(nat),allep(nat))

  ! Initialize atoms positions and give them random velocities
  box=nat
  x=[(real(i,rk),i=begindex-1,endindex-1)]
  
  seed=seed+myid ! so that the random velocities in each process are not equal
  do i=1,niter
     v(i)=lcg(seed) ! filling the velocity array with random numbers
  end do
  v=vsc*v

  ! Remove center of mass velocity

  call mpi_allreduce(sum(v), vsum, 1, mpi_double_precision, &
       &mpi_sum, mpi_comm_world, ierr) ! sum all velocities in all processes
  v=v-vsum/nat ! alters the v array in order to remove the center of mass
 
  n=0

  ! If the user wants calculate initial energy and print initial coords
  if (cout>0) then ! this part of the code is ok - no problems here!
     ! send and receive information about the position of the atoms at the borders
     call neighboor_communication (x(1),x(niter),myid,ntasks,msg)
     do i=1,niter
        call accel(i,ep(i),a(i))
     end do
!     print *, 'Begin gathering.'  ! debugging
     call mpi_gather(x,niter,mpi_double_precision,&
           &allcoord,niter,mpi_double_precision,ROOT,mpi_comm_world,ierr)
     call mpi_gather(ep,niter,mpi_double_precision,&
          &allep,niter,mpi_double_precision,ROOT,mpi_comm_world,ierr)
!     print *, 'End gathering.'  ! debugging
     if (myid==ROOT) then ! only the root print things
!        print *, 'I am root. Begin print' !  debugging
        call printcoords()
!        print *, 'End Print' ! debugging
     end if
  end if

  ! Simulation proper

  time_loop: do n=1,maxt

     v0=v
     ! send and receive information about the position of the atoms at the borders
     call neighboor_communication (x(1),x(niter),myid,ntasks,msg)
     atom_loop1: do i=1,niter
        ! New potential energy and acceleration
        call accel(i,ep(i),a(i))
     end do atom_loop1
     
     atom_loop2: do i=1,niter
        
        ! Leap frog integration algorithm: update position and velocity
        v(i)=v(i)+dt*a(i)
        x(i)=x(i)+dt*v(i)

        ! Check periodic boundary conditions at the limits of the box
        if (x(i)<0.0 ) x(i)=x(i)+box
        if (x(i)>=box) x(i)=x(i)-box
        
        ! Calculate kinetic energy (note: mass=1)
        ek(i)=1.0/2.0*((v0(i)+v(i))/2.0)**2
        
     end do atom_loop2

     ! Calculate and print total potential end kinetic energies
     ! and their sum that should be conserved.
     call mpi_reduce(sum(ep), epsum, 1, mpi_double_precision, &
       &mpi_sum, ROOT, mpi_comm_world, ierr) ! sum all the values of epsum and
                                             ! send it to ROOT
     call mpi_reduce(sum(ek), eksum, 1, mpi_double_precision, &
       &mpi_sum, ROOT, mpi_comm_world, ierr) ! the same as above

     if (eout>0) then
        if (myid==ROOT) then ! only the ROOT prints the values 
           if (mod(n,eout)==0) print '(4g20.10)',dt*n,&
                &epsum+eksum,epsum,eksum
        end if
     end if
     if (cout>0) then
        if (mod(n,cout)==0) then
!           print *, 'mod(n,cout)==0'  ! debugging
!           print *, 'begin gather n>0'  ! debugging 
           call mpi_gather(x,niter,mpi_double_precision,&
                &allcoord,niter,mpi_double_precision,ROOT,mpi_comm_world,ierr)
!           print *, 'Middle'  ! debugging 
           call mpi_gather(ep,niter,mpi_double_precision,&
                &allep,niter,mpi_double_precision,ROOT,mpi_comm_world,ierr)
!           print *, 'end gather n>0.'  ! debugging 
           if (myid==ROOT) then
!              print *, 'I am root'  ! debugging 
              call printcoords()
!              print *, 'End print:', n  ! debugging 
           end if
        end if
     end if
     
  end do time_loop

  t1=MPI_Wtime()

  tused=t1-t0
  if (myid==ROOT) then
     if (time_print) print*, tused
  end if

  call mpi_finalize(ierr)
  
  stop

  contains

    real(rk) function lcg(seed)
      implicit none
      integer(ik),intent(inout) :: seed
      integer(ik) :: a=69069,c=1,m=4294967296_ik
      real(rk) :: rm=4294967296.0
      seed=mod((seed*a+c),m)
      lcg=seed/rm
    end function lcg


    subroutine accel(i,u,a)
      ! Calculate the potential energy u 
      ! and acceleration a of atom i.
      integer,intent(in) :: i
      real(rk),intent(out) :: u,a
      real(rk),parameter :: d=1.0,k1=1.0,k2=0.5
      integer :: j,k
      real(rk) :: dxl,dxr

      ! checks if the atom i
      if (i==1) then ! is the first atom of the array
         dxl = x(1) - msg(1) ! msg(1) is the position of the atom i-1
         dxr = x(2) - x(1)
      else if (i==niter) then ! is the last atom of the array
         dxl = x(niter) - x(niter-1)
         dxr =  msg(2) - x(niter) ! msg(2) is the position of the atom i+1
      else ! if the atom is at the middle of the array
         j=i-1
         k=i+1
         dxl = x(i) - x(j)
         dxr = x(k) - x(i)
      end if
      
      if (dxl<-box/2.0) dxl=dxl+box
      if (dxl>=box/2.0) dxl=dxl-box
      if (dxr<-box/2.0) dxr=dxr+box
      if (dxr>=box/2.0) dxr=dxr-box
      dxl=dxl-d;
      dxr=dxr-d;

      u=(k1*(dxl**2+dxr**2)+k2*(dxl**3+dxr**3))/2.0
      a=-(2.0*k1*(dxl-dxr)+3.0*k2*(dxl**2-dxr**2))

      return
    end subroutine accel

    subroutine printcoords()
      integer :: ia
      real(rk),parameter :: xsc=2.35
      write(10,*) nat
      write(10,'(a,x,i0,x,a,3f14.4)') 'timestep ',n,'boxsize',xsc*box,0.0,0.0
      do ia=1,nat
         write(10,'(a,x,4g20.10)') 'Si',xsc*allcoord(ia),0.0,0.0,allep(ia)
      end do
      return
    end subroutine printcoords
      
  end program md1d
