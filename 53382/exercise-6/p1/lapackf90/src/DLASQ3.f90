      SUBROUTINE DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, &
                         ITER, NDIV, IEEE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     May 17, 2000
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            I0, ITER, N0, NDIV, NFAIL, PP
      DOUBLE PRECISION   DESIG, DMIN, QMAX, SIGMA
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASQ3 checks for deflation, computes a shift (TAU) and calls dqds.
!  In case of failure it changes shifts, and tries again until output
!  is positive.
!
!  Arguments
!  =========
!
!  I0     (input) INTEGER
!         First index.
!
!  N0     (input) INTEGER
!         Last index.
!
!  Z      (input) DOUBLE PRECISION array, dimension ( 4*N )
!         Z holds the qd array.
!
!  PP     (input) INTEGER
!         PP=0 for ping, PP=1 for pong.
!
!  DMIN   (output) DOUBLE PRECISION
!         Minimum value of d.
!
!  SIGMA  (output) DOUBLE PRECISION
!         Sum of shifts used in current segment.
!
!  DESIG  (input/output) DOUBLE PRECISION
!         Lower order part of SIGMA
!
!  QMAX   (input) DOUBLE PRECISION
!         Maximum value of q.
!
!  NFAIL  (output) INTEGER
!         Number of times shift was too big.
!
!  ITER   (output) INTEGER
!         Number of iterations.
!
!  NDIV   (output) INTEGER
!         Number of divisions.
!
!  TTYPE  (output) INTEGER
!         Shift type.
!
!  IEEE   (input) LOGICAL
!         Flag for IEEE or non IEEE arithmetic (passed to DLASQ5).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   CBIAS
      PARAMETER          ( CBIAS = 1.50D0 )
      DOUBLE PRECISION   ZERO, QURTR, HALF, ONE, TWO, HUNDRD
      PARAMETER          ( ZERO = 0.0D0, QURTR = 0.250D0, HALF = 0.5D0, &
                           ONE = 1.0D0, TWO = 2.0D0, HUNDRD = 100.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IPN4, J4, N0IN, NN, TTYPE
      DOUBLE PRECISION   DMIN1, DMIN2, DN, DN1, DN2, EPS, S, SAFMIN, T, &
                         TAU, TEMP, TOL, TOL2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASQ4, DLASQ5, DLASQ6
!     ..
!     .. External Function ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN, SQRT
!     ..
!     .. Save statement ..
      SAVE               TTYPE
      SAVE               DMIN1, DMIN2, DN, DN1, DN2, TAU
!     ..
!     .. Data statement ..
      DATA               TTYPE / 0 /
      DATA               DMIN1 / ZERO /, DMIN2 / ZERO /, DN / ZERO /, &
                         DN1 / ZERO /, DN2 / ZERO /, TAU / ZERO /
!     ..
!     .. Executable Statements ..
!
      N0IN = N0
      EPS = DLAMCH( 'Precision' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      TOL = EPS*HUNDRD
      TOL2 = TOL**2
!
!     Check for deflation.
!
   10 CONTINUE
!
      if ( N0 < I0 ) &
         RETURN
      if ( N0 == I0 ) &
         GO TO 20
      NN = 4*N0 + PP
      if ( N0 == ( I0+1 ) ) &
         GO TO 40
!
!     Check whether E(N0-1) is negligible, 1 eigenvalue.
!
      if ( Z( NN-5 ).GT.TOL2*( SIGMA+Z( NN-3 ) ) .AND. &
          Z( NN-2*PP-4 ).GT.TOL2*Z( NN-7 ) ) &
         GO TO 30
!
   20 CONTINUE
!
      Z( 4*N0-3 ) = Z( 4*N0+PP-3 ) + SIGMA
      N0 = N0 - 1
      GO TO 10
!
!     Check  whether E(N0-2) is negligible, 2 eigenvalues.
!
   30 CONTINUE
!
      if ( Z( NN-9 ).GT.TOL2*SIGMA .AND. &
          Z( NN-2*PP-8 ).GT.TOL2*Z( NN-11 ) ) &
         GO TO 50
!
   40 CONTINUE
!
      if ( Z( NN-3 ).GT.Z( NN-7 ) ) then
         S = Z( NN-3 )
         Z( NN-3 ) = Z( NN-7 )
         Z( NN-7 ) = S
      end if
      if ( Z( NN-5 ).GT.Z( NN-3 )*TOL2 ) then
         T = HALF*( ( Z( NN-7 )-Z( NN-3 ) )+Z( NN-5 ) )
         S = Z( NN-3 )*( Z( NN-5 ) / T )
         if ( S.LE.T ) then
            S = Z( NN-3 )*( Z( NN-5 ) / &
                ( T*( ONE+SQRT( ONE+S / T ) ) ) )
         ELSE
            S = Z( NN-3 )*( Z( NN-5 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
         end if
         T = Z( NN-7 ) + ( S+Z( NN-5 ) )
         Z( NN-3 ) = Z( NN-3 )*( Z( NN-7 ) / T )
         Z( NN-7 ) = T
      end if
      Z( 4*N0-7 ) = Z( NN-7 ) + SIGMA
      Z( 4*N0-3 ) = Z( NN-3 ) + SIGMA
      N0 = N0 - 2
      GO TO 10
!
   50 CONTINUE
!
!     Reverse the qd-array, if warranted.
!
      if ( DMIN.LE.ZERO .OR. N0 < N0IN ) then
         if ( CBIAS*Z( 4*I0+PP-3 ) < Z( 4*N0+PP-3 ) ) then
            IPN4 = 4*( I0+N0 )
            DO 60 J4 = 4*I0, 2*( I0+N0-1 ), 4
               TEMP = Z( J4-3 )
               Z( J4-3 ) = Z( IPN4-J4-3 )
               Z( IPN4-J4-3 ) = TEMP
               TEMP = Z( J4-2 )
               Z( J4-2 ) = Z( IPN4-J4-2 )
               Z( IPN4-J4-2 ) = TEMP
               TEMP = Z( J4-1 )
               Z( J4-1 ) = Z( IPN4-J4-5 )
               Z( IPN4-J4-5 ) = TEMP
               TEMP = Z( J4 )
               Z( J4 ) = Z( IPN4-J4-4 )
               Z( IPN4-J4-4 ) = TEMP
   60       CONTINUE
            if ( N0-I0.LE.4 ) then
               Z( 4*N0+PP-1 ) = Z( 4*I0+PP-1 )
               Z( 4*N0-PP ) = Z( 4*I0-PP )
            end if
            DMIN2 = MIN( DMIN2, Z( 4*N0+PP-1 ) )
            Z( 4*N0+PP-1 ) = MIN( Z( 4*N0+PP-1 ), Z( 4*I0+PP-1 ), &
                                  Z( 4*I0+PP+3 ) )
            Z( 4*N0-PP ) = MIN( Z( 4*N0-PP ), Z( 4*I0-PP ), &
                                Z( 4*I0-PP+4 ) )
            QMAX = MAX( QMAX, Z( 4*I0+PP-3 ), Z( 4*I0+PP+1 ) )
            DMIN = -ZERO
         end if
      end if
!
!  70 CONTINUE
!
      if ( DMIN < ZERO .OR. SAFMIN*QMAX.LT.MIN( Z( 4*N0+PP-1 ), &
          Z( 4*N0+PP-9 ), DMIN2+Z( 4*N0-PP ) ) ) then
!
!        Choose a shift.
!
         CALL DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1, &
                      DN2, TAU, TTYPE )
!
!        Call dqds until DMIN > 0.
!
   80    CONTINUE
!
         CALL DLASQ5( I0, N0, Z, PP, TAU, DMIN, DMIN1, DMIN2, DN, &
                      DN1, DN2, IEEE )
!
         NDIV = NDIV + ( N0-I0+2 )
         ITER = ITER + 1
!
!        Check status.
!
         if ( DMIN.GE.ZERO .AND. DMIN1.GT.ZERO ) then
!
!           Success.
!
            GO TO 100
!
         else if ( DMIN < ZERO .AND. DMIN1.GT.ZERO .AND. &
                  Z( 4*( N0-1 )-PP ) < TOL*( SIGMA+DN1 ) .AND. &
                  ABS( DN ) < TOL*SIGMA ) then
!
!           Convergence hidden by negative DN.
!
            Z( 4*( N0-1 )-PP+2 ) = ZERO
            DMIN = ZERO
            GO TO 100
         else if ( DMIN < ZERO ) then
!
!           TAU too big. Select new TAU and try again.
!
            NFAIL = NFAIL + 1
            if ( TTYPE < -22 ) then
!
!              Failed twice. Play it safe.
!
               TAU = ZERO
            else if ( DMIN1.GT.ZERO ) then
!
!              Late failure. Gives excellent shift.
!
               TAU = ( TAU+DMIN )*( ONE-TWO*EPS )
               TTYPE = TTYPE - 11
            ELSE
!
!              Early failure. Divide by 4.
!
               TAU = QURTR*TAU
               TTYPE = TTYPE - 12
            end if
            GO TO 80
         else if ( DMIN.NE.DMIN ) then
!
!           NaN.
!
            TAU = ZERO
            GO TO 80
         ELSE
!
!           Possible underflow. Play it safe.
!
            GO TO 90
         end if
      end if
!
!     Risk of underflow.
!
   90 CONTINUE
      CALL DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DN1, DN2 )
      NDIV = NDIV + ( N0-I0+2 )
      ITER = ITER + 1
      TAU = ZERO
!
  100 CONTINUE
      if ( TAU < SIGMA ) then
         DESIG = DESIG + TAU
         T = SIGMA + DESIG
         DESIG = DESIG - ( T-SIGMA )
      ELSE
         T = SIGMA + TAU
         DESIG = SIGMA - ( T-TAU ) + DESIG
      end if
      SIGMA = T
!
      RETURN
!
!     End of DLASQ3
!
      END
