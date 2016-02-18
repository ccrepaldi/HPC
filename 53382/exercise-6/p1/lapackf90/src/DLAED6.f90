      SUBROUTINE DLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
!     Courant Institute, NAG Ltd., and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      LOGICAL            ORGATI
      INTEGER            INFO, KNITER
      DOUBLE PRECISION   FINIT, RHO, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( 3 ), Z( 3 )
!     ..
!
!  Purpose
!  =======
!
!  DLAED6 computes the positive or negative root (closest to the origin)
!  of
!                   z(1)        z(2)        z(3)
!  f(x) =   rho + --------- + ---------- + ---------
!                  d(1)-x      d(2)-x      d(3)-x
!
!  It is assumed that
!
!        if ORGATI = .true. the root is between d(2) and d(3);
!        otherwise it is between d(1) and d(2)
!
!  This routine will be called by DLAED4 when necessary. In most cases,
!  the root sought is the smallest in magnitude, though it might not be
!  in some extremely rare situations.
!
!  Arguments
!  =========
!
!  KNITER       (input) INTEGER
!               Refer to DLAED4 for its significance.
!
!  ORGATI       (input) LOGICAL
!               If ORGATI is true, the needed root is between d(2) and
!               d(3); otherwise it is between d(1) and d(2).  See
!               DLAED4 for further details.
!
!  RHO          (input) DOUBLE PRECISION
!               Refer to the equation f(x) above.
!
!  D            (input) DOUBLE PRECISION array, dimension (3)
!               D satisfies d(1) < d(2) < d(3).
!
!  Z            (input) DOUBLE PRECISION array, dimension (3)
!               Each of the elements in z must be positive.
!
!  FINIT        (input) DOUBLE PRECISION
!               The value of f at 0. It is more accurate than the one
!               evaluated inside this routine (if someone wants to do
!               so).
!
!  TAU          (output) DOUBLE PRECISION
!               The root of the equation f(x).
!
!  INFO         (output) INTEGER
!               = 0: successful exit
!               > 0: if INFO = 1, failure to converge
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Ren-Cang Li, Computer Science Division, University of California
!     at Berkeley, USA
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 20 )
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR, EIGHT
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                         THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0 )
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   DSCALE( 3 ), ZSCALE( 3 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FIRST, SCALE
      INTEGER            I, ITER, NITER
      DOUBLE PRECISION   A, B, BASE, C, DDF, DF, EPS, ERRETM, ETA, F, &
                         FC, SCLFAC, SCLINV, SMALL1, SMALL2, SMINV1, &
                         SMINV2, TEMP, TEMP1, TEMP2, TEMP3, TEMP4
!     ..
!     .. Save statement ..
      SAVE               FIRST, SMALL1, SMINV1, SMALL2, SMINV2, EPS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
      NITER = 1
      TAU = ZERO
      if ( KNITER == 2 ) then
         if ( ORGATI ) then
            TEMP = ( D( 3 )-D( 2 ) ) / TWO
            C = RHO + Z( 1 ) / ( ( D( 1 )-D( 2 ) )-TEMP )
            A = C*( D( 2 )+D( 3 ) ) + Z( 2 ) + Z( 3 )
            B = C*D( 2 )*D( 3 ) + Z( 2 )*D( 3 ) + Z( 3 )*D( 2 )
         ELSE
            TEMP = ( D( 1 )-D( 2 ) ) / TWO
            C = RHO + Z( 3 ) / ( ( D( 3 )-D( 2 ) )-TEMP )
            A = C*( D( 1 )+D( 2 ) ) + Z( 1 ) + Z( 2 )
            B = C*D( 1 )*D( 2 ) + Z( 1 )*D( 2 ) + Z( 2 )*D( 1 )
         end if
         TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
         A = A / TEMP
         B = B / TEMP
         C = C / TEMP
         if ( C == ZERO ) then
            TAU = B / A
         else if ( A.LE.ZERO ) then
            TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
         end if
         TEMP = RHO + Z( 1 ) / ( D( 1 )-TAU ) + &
                Z( 2 ) / ( D( 2 )-TAU ) + Z( 3 ) / ( D( 3 )-TAU )
         if ( ABS( FINIT ).LE.ABS( TEMP ) ) &
            TAU = ZERO
      end if
!
!     On first call to routine, get machine parameters for
!     possible scaling to avoid overflow
!
      if ( FIRST ) then
         EPS = DLAMCH( 'Epsilon' )
         BASE = DLAMCH( 'Base' )
         SMALL1 = BASE**( INT( LOG( DLAMCH( 'SafMin' ) ) / LOG( BASE ) / &
                  THREE ) )
         SMINV1 = ONE / SMALL1
         SMALL2 = SMALL1*SMALL1
         SMINV2 = SMINV1*SMINV1
         FIRST = .FALSE.
      end if
!
!     Determine if scaling of inputs necessary to avoid overflow
!     when computing 1/TEMP**3
!
      if ( ORGATI ) then
         TEMP = MIN( ABS( D( 2 )-TAU ), ABS( D( 3 )-TAU ) )
      ELSE
         TEMP = MIN( ABS( D( 1 )-TAU ), ABS( D( 2 )-TAU ) )
      end if
      SCALE = .FALSE.
      if ( TEMP.LE.SMALL1 ) then
         SCALE = .TRUE.
         if ( TEMP.LE.SMALL2 ) then
!
!        Scale up by power of radix nearest 1/SAFMIN**(2/3)
!
            SCLFAC = SMINV2
            SCLINV = SMALL2
         ELSE
!
!        Scale up by power of radix nearest 1/SAFMIN**(1/3)
!
            SCLFAC = SMINV1
            SCLINV = SMALL1
         end if
!
!        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1)
!
         DO 10 I = 1, 3
            DSCALE( I ) = D( I )*SCLFAC
            ZSCALE( I ) = Z( I )*SCLFAC
   10    CONTINUE
         TAU = TAU*SCLFAC
      ELSE
!
!        Copy D and Z to DSCALE and ZSCALE
!
         DO 20 I = 1, 3
            DSCALE( I ) = D( I )
            ZSCALE( I ) = Z( I )
   20    CONTINUE
      end if
!
      FC = ZERO
      DF = ZERO
      DDF = ZERO
      DO 30 I = 1, 3
         TEMP = ONE / ( DSCALE( I )-TAU )
         TEMP1 = ZSCALE( I )*TEMP
         TEMP2 = TEMP1*TEMP
         TEMP3 = TEMP2*TEMP
         FC = FC + TEMP1 / DSCALE( I )
         DF = DF + TEMP2
         DDF = DDF + TEMP3
   30 CONTINUE
      F = FINIT + TAU*FC
!
      if ( ABS( F ).LE.ZERO ) &
         GO TO 60
!
!        Iteration begins
!
!     It is not hard to see that
!
!           1) Iterations will go up monotonically
!              if FINIT < 0;
!
!           2) Iterations will go down monotonically
!              if FINIT > 0.
!
      ITER = NITER + 1
!
      DO 50 NITER = ITER, MAXIT
!
         if ( ORGATI ) then
            TEMP1 = DSCALE( 2 ) - TAU
            TEMP2 = DSCALE( 3 ) - TAU
         ELSE
            TEMP1 = DSCALE( 1 ) - TAU
            TEMP2 = DSCALE( 2 ) - TAU
         end if
         A = ( TEMP1+TEMP2 )*F - TEMP1*TEMP2*DF
         B = TEMP1*TEMP2*F
         C = F - ( TEMP1+TEMP2 )*DF + TEMP1*TEMP2*DDF
         TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
         A = A / TEMP
         B = B / TEMP
         C = C / TEMP
         if ( C == ZERO ) then
            ETA = B / A
         else if ( A.LE.ZERO ) then
            ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
         end if
         if ( F*ETA.GE.ZERO ) then
            ETA = -F / DF
         end if
!
         TEMP = ETA + TAU
         if ( ORGATI ) then
            if ( ETA.GT.ZERO .AND. TEMP.GE.DSCALE( 3 ) ) &
               ETA = ( DSCALE( 3 )-TAU ) / TWO
            if ( ETA < ZERO .AND. TEMP.LE.DSCALE( 2 ) ) &
               ETA = ( DSCALE( 2 )-TAU ) / TWO
         ELSE
            if ( ETA.GT.ZERO .AND. TEMP.GE.DSCALE( 2 ) ) &
               ETA = ( DSCALE( 2 )-TAU ) / TWO
            if ( ETA < ZERO .AND. TEMP.LE.DSCALE( 1 ) ) &
               ETA = ( DSCALE( 1 )-TAU ) / TWO
         end if
         TAU = TAU + ETA
!
         FC = ZERO
         ERRETM = ZERO
         DF = ZERO
         DDF = ZERO
         DO 40 I = 1, 3
            TEMP = ONE / ( DSCALE( I )-TAU )
            TEMP1 = ZSCALE( I )*TEMP
            TEMP2 = TEMP1*TEMP
            TEMP3 = TEMP2*TEMP
            TEMP4 = TEMP1 / DSCALE( I )
            FC = FC + TEMP4
            ERRETM = ERRETM + ABS( TEMP4 )
            DF = DF + TEMP2
            DDF = DDF + TEMP3
   40    CONTINUE
         F = FINIT + TAU*FC
         ERRETM = EIGHT*( ABS( FINIT )+ABS( TAU )*ERRETM ) + &
                  ABS( TAU )*DF
         if ( ABS( F ).LE.EPS*ERRETM ) &
            GO TO 60
   50 CONTINUE
      INFO = 1
   60 CONTINUE
!
!     Undo scaling
!
      if ( SCALE ) &
         TAU = TAU*SCLINV
      RETURN
!
!     End of DLAED6
!
      END
