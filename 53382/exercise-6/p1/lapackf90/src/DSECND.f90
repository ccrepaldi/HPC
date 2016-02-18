      DOUBLE PRECISION FUNCTION DSECND( )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!  Purpose
!  =======
!
!  DSECND returns the user time for a process in seconds.
!  This version gets the time from the system function ETIME.
!
! =====================================================================
!
!     .. Local Scalars ..
      double precision T1
!     ..
!     .. Local Arrays ..
!     REAL               TARRAY( 2 )
!     ..
!     .. External Functions ..
      REAL               ETIME
      EXTERNAL           ETIME
!     ..
!     .. Executable Statements ..
!
      call cpu_time ( t1 )
      DSECND = t1
      RETURN
!
!     End of DSECND
!
      END
