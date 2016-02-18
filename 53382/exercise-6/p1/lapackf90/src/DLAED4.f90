      SUBROUTINE DLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
!     Courant Institute, NAG Ltd., and Rice University
!     December 23, 1999
!
!     .. Scalar Arguments ..
      INTEGER            I, INFO, N
      DOUBLE PRECISION   DLAM, RHO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DELTA( * ), Z( * )
!     ..
!
!  Purpose
!  =======
!
!  This subroutine computes the I-th updated eigenvalue of a symmetric
!  rank-one modification to a diagonal matrix whose elements are
!  given in the array d, and that
!
!             D(i) < D(j)  for  i < j
!
!  and that RHO > 0.  This is arranged by the calling routine, and is
!  no loss in generality.  The rank-one modified system is thus
!
!             diag( D )  +  RHO *  Z * Z_transpose.
!
!  where we assume the Euclidean norm of Z is 1.
!
!  The method consists of approximating the rational functions in the
!  secular equation by simpler interpolating rational functions.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         The length of all arrays.
!
!  I      (input) INTEGER
!         The index of the eigenvalue to be computed.  1 <= I <= N.
!
!  D      (input) DOUBLE PRECISION array, dimension (N)
!         The original eigenvalues.  It is assumed that they are in
!         order, D(I) < D(J)  for I < J.
!
!  Z      (input) DOUBLE PRECISION array, dimension (N)
!         The components of the updating vector.
!
!  DELTA  (output) DOUBLE PRECISION array, dimension (N)
!         If N .ne. 1, DELTA contains (D(j) - lambda_I) in its  j-th
!         component.  If N = 1, then DELTA(1) = 1.  The vector DELTA
!         contains the information necessary to construct the
!         eigenvectors.
!
!  RHO    (input) DOUBLE PRECISION
!         The scalar in the symmetric updating formula.
!
!  DLAM   (output) DOUBLE PRECISION
!         The computed lambda_I, the I-th updated eigenvalue.
!
!  INFO   (output) INTEGER
!         = 0:  successful exit
!         > 0:  if INFO = 1, the updating process failed.
!
!  Internal Parameters
!  ===================
!
!  Logical variable ORGATI (origin-at-i?) is used for distinguishing
!  whether D(i) or D(i+1) is treated as the origin.
!
!            ORGATI = .true.    origin at i
!            ORGATI = .false.   origin at i+1
!
!   Logical variable SWTCH3 (switch-for-3-poles?) is for noting
!   if we are working with THREE poles!
!
!   MAXIT is the maximum number of iterations allowed for each
!   eigenvalue.
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
      PARAMETER          ( MAXIT = 30 )
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR, EIGHT, TEN
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                         THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0, &
                         TEN = 10.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ORGATI, SWTCH, SWTCH3
      INTEGER            II, IIM1, IIP1, IP1, ITER, J, NITER
      DOUBLE PRECISION   A, B, C, DEL, DLTLB, DLTUB, DPHI, DPSI, DW, &
                         EPS, ERRETM, ETA, MIDPT, PHI, PREW, PSI, &
                         RHOINV, TAU, TEMP, TEMP1, W
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   ZZ( 3 )
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAED5, DLAED6
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Since this routine is called in an inner loop, we do no argument
!     checking.
!
!     Quick return for N=1 and 2.
!
      INFO = 0
      if ( N == 1 ) then
!
!         Presumably, I=1 upon entry
!
         DLAM = D( 1 ) + RHO*Z( 1 )*Z( 1 )
         DELTA( 1 ) = ONE
         RETURN
      end if
      if ( N == 2 ) then
         CALL DLAED5( I, D, Z, DELTA, RHO, DLAM )
         RETURN
      end if
!
!     Compute machine epsilon
!
      EPS = DLAMCH( 'Epsilon' )
      RHOINV = ONE / RHO
!
!     The case I = N
!
      if ( I == N ) then
!
!        Initialize some basic variables
!
         II = N - 1
         NITER = 1
!
!        Calculate initial guess
!
         MIDPT = RHO / TWO
!
!        If ||Z||_2 is not one, then TEMP should be set to
!        RHO * ||Z||_2^2 / TWO
!
         DO 10 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
   10    CONTINUE
!
         PSI = ZERO
         DO 20 J = 1, N - 2
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
   20    CONTINUE
!
         C = RHOINV + PSI
         W = C + Z( II )*Z( II ) / DELTA( II ) + &
             Z( N )*Z( N ) / DELTA( N )
!
         if ( W.LE.ZERO ) then
            TEMP = Z( N-1 )*Z( N-1 ) / ( D( N )-D( N-1 )+RHO ) + &
                   Z( N )*Z( N ) / RHO
            if ( C.LE.TEMP ) then
               TAU = RHO
            ELSE
               DEL = D( N ) - D( N-1 )
               A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
               B = Z( N )*Z( N )*DEL
               if ( A < ZERO ) then
                  TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
               ELSE
                  TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
               end if
            end if
!
!           It can be proved that
!               D(N)+RHO/2 <= LAMBDA(N) < D(N)+TAU <= D(N)+RHO
!
            DLTLB = MIDPT
            DLTUB = RHO
         ELSE
            DEL = D( N ) - D( N-1 )
            A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
            B = Z( N )*Z( N )*DEL
            if ( A < ZERO ) then
               TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
            ELSE
               TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
            end if
!
!           It can be proved that
!               D(N) < D(N)+TAU < LAMBDA(N) < D(N)+RHO/2
!
            DLTLB = ZERO
            DLTUB = MIDPT
         end if
!
         DO 30 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - TAU
   30    CONTINUE
!
!        Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 40 J = 1, II
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   40    CONTINUE
         ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV + &
                  ABS( TAU )*( DPSI+DPHI )
!
         W = RHOINV + PHI + PSI
!
!        Test for convergence
!
         if ( ABS( W ).LE.EPS*ERRETM ) then
            DLAM = D( I ) + TAU
            GO TO 250
         end if
!
         if ( W.LE.ZERO ) then
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         end if
!
!        Calculate the new step
!
         NITER = NITER + 1
         C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
         A = ( DELTA( N-1 )+DELTA( N ) )*W - &
             DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
         B = DELTA( N-1 )*DELTA( N )*W
         if ( C < ZERO ) &
            C = ABS( C )
         if ( C == ZERO ) then
!          ETA = B/A
!           ETA = RHO - TAU
            ETA = DLTUB - TAU
         else if ( A.GE.ZERO ) then
            ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
         end if
!
!        Note, eta should be positive if w is negative, and
!        eta should be negative otherwise. However,
!        if for some reason caused by roundoff, eta*w > 0,
!        we simply use one Newton step instead. This way
!        will guarantee eta*w < 0.
!
         if ( W*ETA.GT.ZERO ) &
            ETA = -W / ( DPSI+DPHI )
         TEMP = TAU + ETA
         if ( TEMP.GT.DLTUB .OR. TEMP < DLTLB ) then
            if ( W < ZERO ) then
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            end if
         end if
         DO 50 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
   50    CONTINUE
!
         TAU = TAU + ETA
!
!        Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 60 J = 1, II
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   60    CONTINUE
         ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV + &
                  ABS( TAU )*( DPSI+DPHI )
!
         W = RHOINV + PHI + PSI
!
!        Main loop to update the values of the array   DELTA
!
         ITER = NITER + 1
!
         DO 90 NITER = ITER, MAXIT
!
!           Test for convergence
!
            if ( ABS( W ).LE.EPS*ERRETM ) then
               DLAM = D( I ) + TAU
               GO TO 250
            end if
!
            if ( W.LE.ZERO ) then
               DLTLB = MAX( DLTLB, TAU )
            ELSE
               DLTUB = MIN( DLTUB, TAU )
            end if
!
!           Calculate the new step
!
            C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
            A = ( DELTA( N-1 )+DELTA( N ) )*W - &
                DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
            B = DELTA( N-1 )*DELTA( N )*W
            if ( A.GE.ZERO ) then
               ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
            end if
!
!           Note, eta should be positive if w is negative, and
!           eta should be negative otherwise. However,
!           if for some reason caused by roundoff, eta*w > 0,
!           we simply use one Newton step instead. This way
!           will guarantee eta*w < 0.
!
            if ( W*ETA.GT.ZERO ) &
               ETA = -W / ( DPSI+DPHI )
            TEMP = TAU + ETA
            if ( TEMP.GT.DLTUB .OR. TEMP < DLTLB ) then
               if ( W < ZERO ) then
                  ETA = ( DLTUB-TAU ) / TWO
               ELSE
                  ETA = ( DLTLB-TAU ) / TWO
               end if
            end if
            DO 70 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
   70       CONTINUE
!
            TAU = TAU + ETA
!
!           Evaluate PSI and the derivative DPSI
!
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 80 J = 1, II
               TEMP = Z( J ) / DELTA( J )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
   80       CONTINUE
            ERRETM = ABS( ERRETM )
!
!           Evaluate PHI and the derivative DPHI
!
            TEMP = Z( N ) / DELTA( N )
            PHI = Z( N )*TEMP
            DPHI = TEMP*TEMP
            ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV + &
                     ABS( TAU )*( DPSI+DPHI )
!
            W = RHOINV + PHI + PSI
   90    CONTINUE
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
         INFO = 1
         DLAM = D( I ) + TAU
         GO TO 250
!
!        End for the case I = N
!
      ELSE
!
!        The case for I < N
!
         NITER = 1
         IP1 = I + 1
!
!        Calculate initial guess
!
         DEL = D( IP1 ) - D( I )
         MIDPT = DEL / TWO
         DO 100 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
  100    CONTINUE
!
         PSI = ZERO
         DO 110 J = 1, I - 1
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
  110    CONTINUE
!
         PHI = ZERO
         DO 120 J = N, I + 2, -1
            PHI = PHI + Z( J )*Z( J ) / DELTA( J )
  120    CONTINUE
         C = RHOINV + PSI + PHI
         W = C + Z( I )*Z( I ) / DELTA( I ) + &
             Z( IP1 )*Z( IP1 ) / DELTA( IP1 )
!
         if ( W.GT.ZERO ) then
!
!           d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
!
!           We choose d(i) as origin.
!
            ORGATI = .TRUE.
            A = C*DEL + Z( I )*Z( I ) + Z( IP1 )*Z( IP1 )
            B = Z( I )*Z( I )*DEL
            if ( A.GT.ZERO ) then
               TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            ELSE
               TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            end if
            DLTLB = ZERO
            DLTUB = MIDPT
         ELSE
!
!           (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
!
!           We choose d(i+1) as origin.
!
            ORGATI = .FALSE.
            A = C*DEL - Z( I )*Z( I ) - Z( IP1 )*Z( IP1 )
            B = Z( IP1 )*Z( IP1 )*DEL
            if ( A < ZERO ) then
               TAU = TWO*B / ( A-SQRT( ABS( A*A+FOUR*B*C ) ) )
            ELSE
               TAU = -( A+SQRT( ABS( A*A+FOUR*B*C ) ) ) / ( TWO*C )
            end if
            DLTLB = -MIDPT
            DLTUB = ZERO
         end if
!
         if ( ORGATI ) then
            DO 130 J = 1, N
               DELTA( J ) = ( D( J )-D( I ) ) - TAU
  130       CONTINUE
         ELSE
            DO 140 J = 1, N
               DELTA( J ) = ( D( J )-D( IP1 ) ) - TAU
  140       CONTINUE
         end if
         if ( ORGATI ) then
            II = I
         ELSE
            II = I + 1
         end if
         IIM1 = II - 1
         IIP1 = II + 1
!
!        Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 150 J = 1, IIM1
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  150    CONTINUE
         ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
         DPHI = ZERO
         PHI = ZERO
         DO 160 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  160    CONTINUE
!
         W = RHOINV + PHI + PSI
!
!        W is the value of the secular function with
!        its ii-th element removed.
!
         SWTCH3 = .FALSE.
         if ( ORGATI ) then
            if ( W < ZERO ) &
               SWTCH3 = .TRUE.
         ELSE
            if ( W.GT.ZERO ) &
               SWTCH3 = .TRUE.
         end if
         if ( II == 1 .OR. II.EQ.N ) &
            SWTCH3 = .FALSE.
!
         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = W + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + &
                  THREE*ABS( TEMP ) + ABS( TAU )*DW
!
!        Test for convergence
!
         if ( ABS( W ).LE.EPS*ERRETM ) then
            if ( ORGATI ) then
               DLAM = D( I ) + TAU
            ELSE
               DLAM = D( IP1 ) + TAU
            end if
            GO TO 250
         end if
!
         if ( W.LE.ZERO ) then
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         end if
!
!        Calculate the new step
!
         NITER = NITER + 1
         if ( .NOT.SWTCH3 ) then
            if ( ORGATI ) then
               C = W - DELTA( IP1 )*DW - ( D( I )-D( IP1 ) )* &
                   ( Z( I ) / DELTA( I ) )**2
            ELSE
               C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )* &
                   ( Z( IP1 ) / DELTA( IP1 ) )**2
            end if
            A = ( DELTA( I )+DELTA( IP1 ) )*W - &
                DELTA( I )*DELTA( IP1 )*DW
            B = DELTA( I )*DELTA( IP1 )*W
            if ( C == ZERO ) then
               if ( A == ZERO ) then
                  if ( ORGATI ) then
                     A = Z( I )*Z( I ) + DELTA( IP1 )*DELTA( IP1 )* &
                         ( DPSI+DPHI )
                  ELSE
                     A = Z( IP1 )*Z( IP1 ) + DELTA( I )*DELTA( I )* &
                         ( DPSI+DPHI )
                  end if
               end if
               ETA = B / A
            else if ( A.LE.ZERO ) then
               ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            end if
         ELSE
!
!           Interpolation using THREE most relevant poles
!
            TEMP = RHOINV + PSI + PHI
            if ( ORGATI ) then
               TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) - &
                   ( D( IIM1 )-D( IIP1 ) )*TEMP1
               ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
               ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )* &
                         ( ( DPSI-TEMP1 )+DPHI )
            ELSE
               TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) - &
                   ( D( IIP1 )-D( IIM1 ) )*TEMP1
               ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )* &
                         ( DPSI+( DPHI-TEMP1 ) )
               ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
            end if
            ZZ( 2 ) = Z( II )*Z( II )
            CALL DLAED6( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA, &
                         INFO )
            if ( INFO.NE.0 ) &
               GO TO 250
         end if
!
!        Note, eta should be positive if w is negative, and
!        eta should be negative otherwise. However,
!        if for some reason caused by roundoff, eta*w > 0,
!        we simply use one Newton step instead. This way
!        will guarantee eta*w < 0.
!
         if ( W*ETA.GE.ZERO ) &
            ETA = -W / DW
         TEMP = TAU + ETA
         if ( TEMP.GT.DLTUB .OR. TEMP < DLTLB ) then
            if ( W < ZERO ) then
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            end if
         end if
!
         PREW = W
!
! 170    CONTINUE
         DO 180 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
  180    CONTINUE
!
!        Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 190 J = 1, IIM1
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  190    CONTINUE
         ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
         DPHI = ZERO
         PHI = ZERO
         DO 200 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  200    CONTINUE
!
         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = RHOINV + PHI + PSI + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + &
                  THREE*ABS( TEMP ) + ABS( TAU+ETA )*DW
!
         SWTCH = .FALSE.
         if ( ORGATI ) then
            if ( -W.GT.ABS( PREW ) / TEN ) &
               SWTCH = .TRUE.
         ELSE
            if ( W.GT.ABS( PREW ) / TEN ) &
               SWTCH = .TRUE.
         end if
!
         TAU = TAU + ETA
!
!        Main loop to update the values of the array   DELTA
!
         ITER = NITER + 1
!
         DO 240 NITER = ITER, MAXIT
!
!           Test for convergence
!
            if ( ABS( W ).LE.EPS*ERRETM ) then
               if ( ORGATI ) then
                  DLAM = D( I ) + TAU
               ELSE
                  DLAM = D( IP1 ) + TAU
               end if
               GO TO 250
            end if
!
            if ( W.LE.ZERO ) then
               DLTLB = MAX( DLTLB, TAU )
            ELSE
               DLTUB = MIN( DLTUB, TAU )
            end if
!
!           Calculate the new step
!
            if ( .NOT.SWTCH3 ) then
               if ( .NOT.SWTCH ) then
                  if ( ORGATI ) then
                     C = W - DELTA( IP1 )*DW - &
                         ( D( I )-D( IP1 ) )*( Z( I ) / DELTA( I ) )**2
                  ELSE
                     C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )* &
                         ( Z( IP1 ) / DELTA( IP1 ) )**2
                  end if
               ELSE
                  TEMP = Z( II ) / DELTA( II )
                  if ( ORGATI ) then
                     DPSI = DPSI + TEMP*TEMP
                  ELSE
                     DPHI = DPHI + TEMP*TEMP
                  end if
                  C = W - DELTA( I )*DPSI - DELTA( IP1 )*DPHI
               end if
               A = ( DELTA( I )+DELTA( IP1 ) )*W - &
                   DELTA( I )*DELTA( IP1 )*DW
               B = DELTA( I )*DELTA( IP1 )*W
               if ( C == ZERO ) then
                  if ( A == ZERO ) then
                     if ( .NOT.SWTCH ) then
                        if ( ORGATI ) then
                           A = Z( I )*Z( I ) + DELTA( IP1 )* &
                               DELTA( IP1 )*( DPSI+DPHI )
                        ELSE
                           A = Z( IP1 )*Z( IP1 ) + &
                               DELTA( I )*DELTA( I )*( DPSI+DPHI )
                        end if
                     ELSE
                        A = DELTA( I )*DELTA( I )*DPSI + &
                            DELTA( IP1 )*DELTA( IP1 )*DPHI
                     end if
                  end if
                  ETA = B / A
               else if ( A.LE.ZERO ) then
                  ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
               ELSE
                  ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
               end if
            ELSE
!
!              Interpolation using THREE most relevant poles
!
               TEMP = RHOINV + PSI + PHI
               if ( SWTCH ) then
                  C = TEMP - DELTA( IIM1 )*DPSI - DELTA( IIP1 )*DPHI
                  ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*DPSI
                  ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*DPHI
               ELSE
                  if ( ORGATI ) then
                     TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) - &
                         ( D( IIM1 )-D( IIP1 ) )*TEMP1
                     ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
                     ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )* &
                               ( ( DPSI-TEMP1 )+DPHI )
                  ELSE
                     TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) - &
                         ( D( IIP1 )-D( IIM1 ) )*TEMP1
                     ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )* &
                               ( DPSI+( DPHI-TEMP1 ) )
                     ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
                  end if
               end if
               CALL DLAED6( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA, &
                            INFO )
               if ( INFO.NE.0 ) &
                  GO TO 250
            end if
!
!           Note, eta should be positive if w is negative, and
!           eta should be negative otherwise. However,
!           if for some reason caused by roundoff, eta*w > 0,
!           we simply use one Newton step instead. This way
!           will guarantee eta*w < 0.
!
            if ( W*ETA.GE.ZERO ) &
               ETA = -W / DW
            TEMP = TAU + ETA
            if ( TEMP.GT.DLTUB .OR. TEMP < DLTLB ) then
               if ( W < ZERO ) then
                  ETA = ( DLTUB-TAU ) / TWO
               ELSE
                  ETA = ( DLTLB-TAU ) / TWO
               end if
            end if
!
            DO 210 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
  210       CONTINUE
!
            TAU = TAU + ETA
            PREW = W
!
!           Evaluate PSI and the derivative DPSI
!
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 220 J = 1, IIM1
               TEMP = Z( J ) / DELTA( J )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
  220       CONTINUE
            ERRETM = ABS( ERRETM )
!
!           Evaluate PHI and the derivative DPHI
!
            DPHI = ZERO
            PHI = ZERO
            DO 230 J = N, IIP1, -1
               TEMP = Z( J ) / DELTA( J )
               PHI = PHI + Z( J )*TEMP
               DPHI = DPHI + TEMP*TEMP
               ERRETM = ERRETM + PHI
  230       CONTINUE
!
            TEMP = Z( II ) / DELTA( II )
            DW = DPSI + DPHI + TEMP*TEMP
            TEMP = Z( II )*TEMP
            W = RHOINV + PHI + PSI + TEMP
            ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + &
                     THREE*ABS( TEMP ) + ABS( TAU )*DW
            if ( W*PREW.GT.ZERO .AND. ABS( W ).GT.ABS( PREW ) / TEN ) &
               SWTCH = .NOT.SWTCH
!
  240    CONTINUE
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
         INFO = 1
         if ( ORGATI ) then
            DLAM = D( I ) + TAU
         ELSE
            DLAM = D( IP1 ) + TAU
         end if
!
      end if
!
  250 CONTINUE
!
      RETURN
!
!     End of DLAED4
!
      END
