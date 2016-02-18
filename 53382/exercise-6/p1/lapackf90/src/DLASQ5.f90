      SUBROUTINE DLASQ5( I0, N0, Z, PP, TAU, DMIN, DMIN1, DMIN2, DN, &
                         DNM1, DNM2, IEEE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     May 17, 2000
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            I0, N0, PP
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASQ5 computes one dqds transform in ping-pong form, one
!  version for IEEE machines another for non IEEE machines.
!
!  Arguments
!  =========
!
!  I0    (input) INTEGER
!        First index.
!
!  N0    (input) INTEGER
!        Last index.
!
!  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )
!        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
!        an extra argument.
!
!  PP    (input) INTEGER
!        PP=0 for ping, PP=1 for pong.
!
!  TAU   (input) DOUBLE PRECISION
!        This is the shift.
!
!  DMIN  (output) DOUBLE PRECISION
!        Minimum value of d.
!
!  DMIN1 (output) DOUBLE PRECISION
!        Minimum value of d, excluding D( N0 ).
!
!  DMIN2 (output) DOUBLE PRECISION
!        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!
!  DN    (output) DOUBLE PRECISION
!        d(N0), the last value of d.
!
!  DNM1  (output) DOUBLE PRECISION
!        d(N0-1).
!
!  DNM2  (output) DOUBLE PRECISION
!        d(N0-2).
!
!  IEEE  (input) LOGICAL
!        Flag for IEEE or non IEEE arithmetic.
!
!  =====================================================================
!
!     .. Parameter ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J4, J4P2
      DOUBLE PRECISION   D, EMIN, TEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      if ( ( N0-I0-1 ).LE.0 ) &
         RETURN
!
      J4 = 4*I0 + PP - 3
      EMIN = Z( J4+4 )
      D = Z( J4 ) - TAU
      DMIN = D
      DMIN1 = -Z( J4 )
!
      if ( IEEE ) then
!
!        Code for IEEE arithmetic.
!
         if ( PP == 0 ) then
            DO 10 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-2 ) = D + Z( J4-1 )
               TEMP = Z( J4+1 ) / Z( J4-2 )
               D = D*TEMP - TAU
               DMIN = MIN( DMIN, D )
               Z( J4 ) = Z( J4-1 )*TEMP
               EMIN = MIN( Z( J4 ), EMIN )
   10       CONTINUE
         ELSE
            DO 20 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-3 ) = D + Z( J4 )
               TEMP = Z( J4+2 ) / Z( J4-3 )
               D = D*TEMP - TAU
               DMIN = MIN( DMIN, D )
               Z( J4-1 ) = Z( J4 )*TEMP
               EMIN = MIN( Z( J4-1 ), EMIN )
   20       CONTINUE
         end if
!
!        Unroll last two steps.
!
         DNM2 = D
         DMIN2 = DMIN
         J4 = 4*( N0-2 ) - PP
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM2 + Z( J4P2 )
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
         DMIN = MIN( DMIN, DNM1 )
!
         DMIN1 = DMIN
         J4 = J4 + 4
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM1 + Z( J4P2 )
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
         DMIN = MIN( DMIN, DN )
!
      ELSE
!
!        Code for non IEEE arithmetic.
!
         if ( PP == 0 ) then
            DO 30 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-2 ) = D + Z( J4-1 )
               if ( D < ZERO ) then
                  RETURN
               ELSE
                  Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
                  D = Z( J4+1 )*( D / Z( J4-2 ) ) - TAU
               end if
               DMIN = MIN( DMIN, D )
               EMIN = MIN( EMIN, Z( J4 ) )
   30       CONTINUE
         ELSE
            DO 40 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-3 ) = D + Z( J4 )
               if ( D < ZERO ) then
                  RETURN
               ELSE
                  Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
                  D = Z( J4+2 )*( D / Z( J4-3 ) ) - TAU
               end if
               DMIN = MIN( DMIN, D )
               EMIN = MIN( EMIN, Z( J4-1 ) )
   40       CONTINUE
         end if
!
!        Unroll last two steps.
!
         DNM2 = D
         DMIN2 = DMIN
         J4 = 4*( N0-2 ) - PP
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM2 + Z( J4P2 )
         if ( DNM2 < ZERO ) then
            RETURN
         ELSE
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
         end if
         DMIN = MIN( DMIN, DNM1 )
!
         DMIN1 = DMIN
         J4 = J4 + 4
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM1 + Z( J4P2 )
         if ( DNM1 < ZERO ) then
            RETURN
         ELSE
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
         end if
         DMIN = MIN( DMIN, DN )
!
      end if
!
      Z( J4+2 ) = DN
      Z( 4*N0-PP ) = EMIN
      RETURN
!
!     End of DLASQ5
!
      END
