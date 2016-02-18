      SUBROUTINE DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, &
                         DN1, DN2, TAU, TTYPE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1999
!
!     .. Scalar Arguments ..
      INTEGER            I0, N0, N0IN, PP, TTYPE
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DN1, DN2, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASQ4 computes an approximation TAU to the smallest eigenvalue
!  using values of d from the previous transform.
!
!  I0    (input) INTEGER
!        First index.
!
!  N0    (input) INTEGER
!        Last index.
!
!  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )
!        Z holds the qd array.
!
!  PP    (input) INTEGER
!        PP=0 for ping, PP=1 for pong.
!
!  NOIN  (input) INTEGER
!        The value of N0 at start of EIGTEST.
!
!  DMIN  (input) DOUBLE PRECISION
!        Minimum value of d.
!
!  DMIN1 (input) DOUBLE PRECISION
!        Minimum value of d, excluding D( N0 ).
!
!  DMIN2 (input) DOUBLE PRECISION
!        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!
!  DN    (input) DOUBLE PRECISION
!        d(N)
!
!  DN1   (input) DOUBLE PRECISION
!        d(N-1)
!
!  DN2   (input) DOUBLE PRECISION
!        d(N-2)
!
!  TAU   (output) DOUBLE PRECISION
!        This is the shift.
!
!  TTYPE (output) INTEGER
!        Shift type.
!
!  Further Details
!  ===============
!  CNST1 = 9/16
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   CNST1, CNST2, CNST3
      PARAMETER          ( CNST1 = 0.5630D0, CNST2 = 1.010D0, &
                         CNST3 = 1.050D0 )
      DOUBLE PRECISION   QURTR, THIRD, HALF, ZERO, ONE, TWO, HUNDRD
      PARAMETER          ( QURTR = 0.250D0, THIRD = 0.3330D0, &
                         HALF = 0.50D0, ZERO = 0.0D0, ONE = 1.0D0, &
                         TWO = 2.0D0, HUNDRD = 100.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I4, NN, NP
      DOUBLE PRECISION   A2, B1, B2, G, GAM, GAP1, GAP2, S
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
!     ..
!     .. Save statement ..
      SAVE               G
!     ..
!     .. Data statement ..
      DATA               G / ZERO /
!     ..
!     .. Executable Statements ..
!
!     A negative DMIN forces the shift to take that absolute value
!     TTYPE records the type of shift.
!
      if ( DMIN.LE.ZERO ) then
         TAU = -DMIN
         TTYPE = -1
         RETURN
      end if
!
      NN = 4*N0 + PP
      if ( N0IN == N0 ) then
!
!        No eigenvalues deflated.
!
         if ( DMIN == DN .OR. DMIN.EQ.DN1 ) then
!
            B1 = SQRT( Z( NN-3 ) )*SQRT( Z( NN-5 ) )
            B2 = SQRT( Z( NN-7 ) )*SQRT( Z( NN-9 ) )
            A2 = Z( NN-7 ) + Z( NN-5 )
!
!           Cases 2 and 3.
!
            if ( DMIN == DN .AND. DMIN1.EQ.DN1 ) then
               GAP2 = DMIN2 - A2 - DMIN2*QURTR
               if ( GAP2.GT.ZERO .AND. GAP2.GT.B2 ) then
                  GAP1 = A2 - DN - ( B2 / GAP2 )*B2
               ELSE
                  GAP1 = A2 - DN - ( B1+B2 )
               end if
               if ( GAP1.GT.ZERO .AND. GAP1.GT.B1 ) then
                  S = MAX( DN-( B1 / GAP1 )*B1, HALF*DMIN )
                  TTYPE = -2
               ELSE
                  S = ZERO
                  if ( DN.GT.B1 ) &
                     S = DN - B1
                  if ( A2.GT.( B1+B2 ) ) &
                     S = MIN( S, A2-( B1+B2 ) )
                  S = MAX( S, THIRD*DMIN )
                  TTYPE = -3
               end if
            ELSE
!
!              Case 4.
!
               TTYPE = -4
               S = QURTR*DMIN
               if ( DMIN == DN ) then
                  GAM = DN
                  A2 = ZERO
                  if ( Z( NN-5 ) .GT. Z( NN-7 ) ) &
                     RETURN
                  B2 = Z( NN-5 ) / Z( NN-7 )
                  NP = NN - 9
               ELSE
                  NP = NN - 2*PP
                  B2 = Z( NP-2 )
                  GAM = DN1
                  if ( Z( NP-4 ) .GT. Z( NP-2 ) ) &
                     RETURN
                  A2 = Z( NP-4 ) / Z( NP-2 )
                  if ( Z( NN-9 ) .GT. Z( NN-11 ) ) &
                     RETURN
                  B2 = Z( NN-9 ) / Z( NN-11 )
                  NP = NN - 13
               end if
!
!              Approximate contribution to norm squared from I < NN-1.
!
               A2 = A2 + B2
               DO 10 I4 = NP, 4*I0 - 1 + PP, -4
                  if ( B2 == ZERO ) &
                     GO TO 20
                  B1 = B2
                  if ( Z( I4 ) .GT. Z( I4-2 ) ) &
                     RETURN
                  B2 = B2*( Z( I4 ) / Z( I4-2 ) )
                  A2 = A2 + B2
                  if ( HUNDRD*MAX( B2, B1 ) < A2 .OR. CNST1.LT.A2 ) &
                     GO TO 20
   10          CONTINUE
   20          CONTINUE
               A2 = CNST3*A2
!
!              Rayleigh quotient residual bound.
!
               if ( A2 < CNST1 ) &
                  S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )
            end if
         else if ( DMIN == DN2 ) then
!
!           Case 5.
!
            TTYPE = -5
            S = QURTR*DMIN
!
!           Compute contribution to norm squared from I > NN-2.
!
            NP = NN - 2*PP
            B1 = Z( NP-2 )
            B2 = Z( NP-6 )
            GAM = DN2
            if ( Z( NP-8 ).GT.B2 .OR. Z( NP-4 ).GT.B1 ) &
               RETURN
            A2 = ( Z( NP-8 ) / B2 )*( ONE+Z( NP-4 ) / B1 )
!
!           Approximate contribution to norm squared from I < NN-2.
!
            if ( N0-I0.GT.2 ) then
               B2 = Z( NN-13 ) / Z( NN-15 )
               A2 = A2 + B2
               DO 30 I4 = NN - 17, 4*I0 - 1 + PP, -4
                  if ( B2 == ZERO ) &
                     GO TO 40
                  B1 = B2
                  if ( Z( I4 ) .GT. Z( I4-2 ) ) &
                     RETURN
                  B2 = B2*( Z( I4 ) / Z( I4-2 ) )
                  A2 = A2 + B2
                  if ( HUNDRD*MAX( B2, B1 ) < A2 .OR. CNST1.LT.A2 ) &
                     GO TO 40
   30          CONTINUE
   40          CONTINUE
               A2 = CNST3*A2
            end if
!
            if ( A2 < CNST1 ) &
               S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )
         ELSE
!
!           Case 6, no information to guide us.
!
            if ( TTYPE == -6 ) then
               G = G + THIRD*( ONE-G )
            else if ( TTYPE == -18 ) then
               G = QURTR*THIRD
            ELSE
               G = QURTR
            end if
            S = G*DMIN
            TTYPE = -6
         end if
!
      else if ( N0IN == ( N0+1 ) ) then
!
!        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.
!
         if ( DMIN1 == DN1 .AND. DMIN2.EQ.DN2 ) then
!
!           Cases 7 and 8.
!
            TTYPE = -7
            S = THIRD*DMIN1
            if ( Z( NN-5 ).GT.Z( NN-7 ) ) &
               RETURN
            B1 = Z( NN-5 ) / Z( NN-7 )
            B2 = B1
            if ( B2 == ZERO ) &
               GO TO 60
            DO 50 I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
               A2 = B1
               if ( Z( I4 ).GT.Z( I4-2 ) ) &
                  RETURN
               B1 = B1*( Z( I4 ) / Z( I4-2 ) )
               B2 = B2 + B1
               if ( HUNDRD*MAX( B1, A2 ) < B2 ) &
                  GO TO 60
   50       CONTINUE
   60       CONTINUE
            B2 = SQRT( CNST3*B2 )
            A2 = DMIN1 / ( ONE+B2**2 )
            GAP2 = HALF*DMIN2 - A2
            if ( GAP2.GT.ZERO .AND. GAP2.GT.B2*A2 ) then
               S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
            ELSE
               S = MAX( S, A2*( ONE-CNST2*B2 ) )
               TTYPE = -8
            end if
         ELSE
!
!           Case 9.
!
            S = QURTR*DMIN1
            if ( DMIN1 == DN1 ) &
               S = HALF*DMIN1
            TTYPE = -9
         end if
!
      else if ( N0IN == ( N0+2 ) ) then
!
!        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.
!
!        Cases 10 and 11.
!
         if ( DMIN2 == DN2 .AND. TWO*Z( NN-5 ) < Z( NN-7 ) ) then
            TTYPE = -10
            S = THIRD*DMIN2
            if ( Z( NN-5 ).GT.Z( NN-7 ) ) &
               RETURN
            B1 = Z( NN-5 ) / Z( NN-7 )
            B2 = B1
            if ( B2 == ZERO ) &
               GO TO 80
            DO 70 I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
               if ( Z( I4 ).GT.Z( I4-2 ) ) &
                  RETURN
               B1 = B1*( Z( I4 ) / Z( I4-2 ) )
               B2 = B2 + B1
               if ( HUNDRD*B1 < B2 ) &
                  GO TO 80
   70       CONTINUE
   80       CONTINUE
            B2 = SQRT( CNST3*B2 )
            A2 = DMIN2 / ( ONE+B2**2 )
            GAP2 = Z( NN-7 ) + Z( NN-9 ) - &
                   SQRT( Z( NN-11 ) )*SQRT( Z( NN-9 ) ) - A2
            if ( GAP2.GT.ZERO .AND. GAP2.GT.B2*A2 ) then
               S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
            ELSE
               S = MAX( S, A2*( ONE-CNST2*B2 ) )
            end if
         ELSE
            S = QURTR*DMIN2
            TTYPE = -11
         end if
      else if ( N0IN.GT.( N0+2 ) ) then
!
!        Case 12, more than two eigenvalues deflated. No information.
!
         S = ZERO
         TTYPE = -12
      end if
!
      TAU = S
      RETURN
!
!     End of DLASQ4
!
      END
