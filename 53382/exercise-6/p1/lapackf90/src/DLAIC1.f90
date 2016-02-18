      SUBROUTINE DLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            J, JOB
      DOUBLE PRECISION   C, GAMMA, S, SEST, SESTPR
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   W( J ), X( J )
!     ..
!
!  Purpose
!  =======
!
!  DLAIC1 applies one step of incremental condition estimation in
!  its simplest version:
!
!  Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j
!  lower triangular matrix L, such that
!           twonorm(L*x) = sest
!  Then DLAIC1 computes sestpr, s, c such that
!  the vector
!                  [ s*x ]
!           xhat = [  c  ]
!  is an approximate singular vector of
!                  [ L     0  ]
!           Lhat = [ w' gamma ]
!  in the sense that
!           twonorm(Lhat*xhat) = sestpr.
!
!  Depending on JOB, an estimate for the largest or smallest singular
!  value is computed.
!
!  Note that [s c]' and sestpr**2 is an eigenpair of the system
!
!      diag(sest*sest, 0) + [alpha  gamma] * [ alpha ]
!                                            [ gamma ]
!
!  where  alpha =  x'*w.
!
!  Arguments
!  =========
!
!  JOB     (input) INTEGER
!          = 1: an estimate for the largest singular value is computed.
!          = 2: an estimate for the smallest singular value is computed.
!
!  J       (input) INTEGER
!          Length of X and W
!
!  X       (input) DOUBLE PRECISION array, dimension (J)
!          The j-vector x.
!
!  SEST    (input) DOUBLE PRECISION
!          Estimated singular value of j by j matrix L
!
!  W       (input) DOUBLE PRECISION array, dimension (J)
!          The j-vector w.
!
!  GAMMA   (input) DOUBLE PRECISION
!          The diagonal element gamma.
!
!  SEDTPR  (output) DOUBLE PRECISION
!          Estimated singular value of (j+1) by (j+1) matrix Lhat.
!
!  S       (output) DOUBLE PRECISION
!          Sine needed in forming xhat.
!
!  C       (output) DOUBLE PRECISION
!          Cosine needed in forming xhat.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      DOUBLE PRECISION   HALF, FOUR
      PARAMETER          ( HALF = 0.5D0, FOUR = 4.0D0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   ABSALP, ABSEST, ABSGAM, ALPHA, B, COSINE, EPS, &
                         NORMA, S1, S2, SINE, T, TEST, TMP, ZETA1, ZETA2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DDOT, DLAMCH
      EXTERNAL           DDOT, DLAMCH
!     ..
!     .. Executable Statements ..
!
      EPS = DLAMCH( 'Epsilon' )
      ALPHA = DDOT( J, X, 1, W, 1 )
!
      ABSALP = ABS( ALPHA )
      ABSGAM = ABS( GAMMA )
      ABSEST = ABS( SEST )
!
      if ( JOB == 1 ) then
!
!        Estimating largest singular value
!
!        special cases
!
         if ( SEST == ZERO ) then
            S1 = MAX( ABSGAM, ABSALP )
            if ( S1 == ZERO ) then
               S = ZERO
               C = ONE
               SESTPR = ZERO
            ELSE
               S = ALPHA / S1
               C = GAMMA / S1
               TMP = SQRT( S*S+C*C )
               S = S / TMP
               C = C / TMP
               SESTPR = S1*TMP
            end if
            RETURN
         else if ( ABSGAM.LE.EPS*ABSEST ) then
            S = ONE
            C = ZERO
            TMP = MAX( ABSEST, ABSALP )
            S1 = ABSEST / TMP
            S2 = ABSALP / TMP
            SESTPR = TMP*SQRT( S1*S1+S2*S2 )
            RETURN
         else if ( ABSALP.LE.EPS*ABSEST ) then
            S1 = ABSGAM
            S2 = ABSEST
            if ( S1.LE.S2 ) then
               S = ONE
               C = ZERO
               SESTPR = S2
            ELSE
               S = ZERO
               C = ONE
               SESTPR = S1
            end if
            RETURN
         else if ( ABSEST.LE.EPS*ABSALP .OR. ABSEST.LE.EPS*ABSGAM ) then
            S1 = ABSGAM
            S2 = ABSALP
            if ( S1.LE.S2 ) then
               TMP = S1 / S2
               S = SQRT( ONE+TMP*TMP )
               SESTPR = S2*S
               C = ( GAMMA / S2 ) / S
               S = SIGN( ONE, ALPHA ) / S
            ELSE
               TMP = S2 / S1
               C = SQRT( ONE+TMP*TMP )
               SESTPR = S1*C
               S = ( ALPHA / S1 ) / C
               C = SIGN( ONE, GAMMA ) / C
            end if
            RETURN
         ELSE
!
!           normal case
!
            ZETA1 = ALPHA / ABSEST
            ZETA2 = GAMMA / ABSEST
!
            B = ( ONE-ZETA1*ZETA1-ZETA2*ZETA2 )*HALF
            C = ZETA1*ZETA1
            if ( B.GT.ZERO ) then
               T = C / ( B+SQRT( B*B+C ) )
            ELSE
               T = SQRT( B*B+C ) - B
            end if
!
            SINE = -ZETA1 / T
            COSINE = -ZETA2 / ( ONE+T )
            TMP = SQRT( SINE*SINE+COSINE*COSINE )
            S = SINE / TMP
            C = COSINE / TMP
            SESTPR = SQRT( T+ONE )*ABSEST
            RETURN
         end if
!
      else if ( JOB == 2 ) then
!
!        Estimating smallest singular value
!
!        special cases
!
         if ( SEST == ZERO ) then
            SESTPR = ZERO
            if ( MAX( ABSGAM, ABSALP ) == ZERO ) then
               SINE = ONE
               COSINE = ZERO
            ELSE
               SINE = -GAMMA
               COSINE = ALPHA
            end if
            S1 = MAX( ABS( SINE ), ABS( COSINE ) )
            S = SINE / S1
            C = COSINE / S1
            TMP = SQRT( S*S+C*C )
            S = S / TMP
            C = C / TMP
            RETURN
         else if ( ABSGAM.LE.EPS*ABSEST ) then
            S = ZERO
            C = ONE
            SESTPR = ABSGAM
            RETURN
         else if ( ABSALP.LE.EPS*ABSEST ) then
            S1 = ABSGAM
            S2 = ABSEST
            if ( S1.LE.S2 ) then
               S = ZERO
               C = ONE
               SESTPR = S1
            ELSE
               S = ONE
               C = ZERO
               SESTPR = S2
            end if
            RETURN
         else if ( ABSEST.LE.EPS*ABSALP .OR. ABSEST.LE.EPS*ABSGAM ) then
            S1 = ABSGAM
            S2 = ABSALP
            if ( S1.LE.S2 ) then
               TMP = S1 / S2
               C = SQRT( ONE+TMP*TMP )
               SESTPR = ABSEST*( TMP / C )
               S = -( GAMMA / S2 ) / C
               C = SIGN( ONE, ALPHA ) / C
            ELSE
               TMP = S2 / S1
               S = SQRT( ONE+TMP*TMP )
               SESTPR = ABSEST / S
               C = ( ALPHA / S1 ) / S
               S = -SIGN( ONE, GAMMA ) / S
            end if
            RETURN
         ELSE
!
!           normal case
!
            ZETA1 = ALPHA / ABSEST
            ZETA2 = GAMMA / ABSEST
!
            NORMA = MAX( ONE+ZETA1*ZETA1+ABS( ZETA1*ZETA2 ), &
                    ABS( ZETA1*ZETA2 )+ZETA2*ZETA2 )
!
!           See if root is closer to zero or to ONE
!
            TEST = ONE + TWO*( ZETA1-ZETA2 )*( ZETA1+ZETA2 )
            if ( TEST.GE.ZERO ) then
!
!              root is close to zero, compute directly
!
               B = ( ZETA1*ZETA1+ZETA2*ZETA2+ONE )*HALF
               C = ZETA2*ZETA2
               T = C / ( B+SQRT( ABS( B*B-C ) ) )
               SINE = ZETA1 / ( ONE-T )
               COSINE = -ZETA2 / T
               SESTPR = SQRT( T+FOUR*EPS*EPS*NORMA )*ABSEST
            ELSE
!
!              root is closer to ONE, shift by that amount
!
               B = ( ZETA2*ZETA2+ZETA1*ZETA1-ONE )*HALF
               C = ZETA1*ZETA1
               if ( B.GE.ZERO ) then
                  T = -C / ( B+SQRT( B*B+C ) )
               ELSE
                  T = B - SQRT( B*B+C )
               end if
               SINE = -ZETA1 / T
               COSINE = -ZETA2 / ( ONE+T )
               SESTPR = SQRT( ONE+T+FOUR*EPS*EPS*NORMA )*ABSEST
            end if
            TMP = SQRT( SINE*SINE+COSINE*COSINE )
            S = SINE / TMP
            C = COSINE / TMP
            RETURN
!
         end if
      end if
      RETURN
!
!     End of DLAIC1
!
      END
