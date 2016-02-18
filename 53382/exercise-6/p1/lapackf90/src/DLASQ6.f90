      SUBROUTINE DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, &
                         DNM1, DNM2 )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1999
!
!     .. Scalar Arguments ..
      INTEGER            I0, N0, PP
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASQ6 computes one dqd (shift equal to zero) transform in
!  ping-pong form, with protection against underflow and overflow.
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
!  =====================================================================
!
!     .. Parameter ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J4, J4P2
      DOUBLE PRECISION   D, EMIN, SAFMIN, TEMP
!     ..
!     .. External Function ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      if ( ( N0-I0-1 ).LE.0 ) &
         RETURN
!
      SAFMIN = DLAMCH( 'Safe minimum' )
      J4 = 4*I0 + PP - 3
      EMIN = Z( J4+4 )
      D = Z( J4 )
      DMIN = D
!
      if ( PP == 0 ) then
         DO 10 J4 = 4*I0, 4*( N0-3 ), 4
            Z( J4-2 ) = D + Z( J4-1 )
            if ( Z( J4-2 ) == ZERO ) then
               Z( J4 ) = ZERO
               D = Z( J4+1 )
               DMIN = D
               EMIN = ZERO
            else if ( SAFMIN*Z( J4+1 ) < Z( J4-2 ) .AND. &
                     SAFMIN*Z( J4-2 ) < Z( J4+1 ) ) then
               TEMP = Z( J4+1 ) / Z( J4-2 )
               Z( J4 ) = Z( J4-1 )*TEMP
               D = D*TEMP
            ELSE
               Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
               D = Z( J4+1 )*( D / Z( J4-2 ) )
            end if
            DMIN = MIN( DMIN, D )
            EMIN = MIN( EMIN, Z( J4 ) )
   10    CONTINUE
      ELSE
         DO 20 J4 = 4*I0, 4*( N0-3 ), 4
            Z( J4-3 ) = D + Z( J4 )
            if ( Z( J4-3 ) == ZERO ) then
               Z( J4-1 ) = ZERO
               D = Z( J4+2 )
               DMIN = D
               EMIN = ZERO
            else if ( SAFMIN*Z( J4+2 ) < Z( J4-3 ) .AND. &
                     SAFMIN*Z( J4-3 ) < Z( J4+2 ) ) then
               TEMP = Z( J4+2 ) / Z( J4-3 )
               Z( J4-1 ) = Z( J4 )*TEMP
               D = D*TEMP
            ELSE
               Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
               D = Z( J4+2 )*( D / Z( J4-3 ) )
            end if
            DMIN = MIN( DMIN, D )
            EMIN = MIN( EMIN, Z( J4-1 ) )
   20    CONTINUE
      end if
!
!     Unroll last two steps.
!
      DNM2 = D
      DMIN2 = DMIN
      J4 = 4*( N0-2 ) - PP
      J4P2 = J4 + 2*PP - 1
      Z( J4-2 ) = DNM2 + Z( J4P2 )
      if ( Z( J4-2 ) == ZERO ) then
         Z( J4 ) = ZERO
         DNM1 = Z( J4P2+2 )
         DMIN = DNM1
         EMIN = ZERO
      else if ( SAFMIN*Z( J4P2+2 ) < Z( J4-2 ) .AND. &
               SAFMIN*Z( J4-2 ) < Z( J4P2+2 ) ) then
         TEMP = Z( J4P2+2 ) / Z( J4-2 )
         Z( J4 ) = Z( J4P2 )*TEMP
         DNM1 = DNM2*TEMP
      ELSE
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) )
      end if
      DMIN = MIN( DMIN, DNM1 )
!
      DMIN1 = DMIN
      J4 = J4 + 4
      J4P2 = J4 + 2*PP - 1
      Z( J4-2 ) = DNM1 + Z( J4P2 )
      if ( Z( J4-2 ) == ZERO ) then
         Z( J4 ) = ZERO
         DN = Z( J4P2+2 )
         DMIN = DN
         EMIN = ZERO
      else if ( SAFMIN*Z( J4P2+2 ) < Z( J4-2 ) .AND. &
               SAFMIN*Z( J4-2 ) < Z( J4P2+2 ) ) then
         TEMP = Z( J4P2+2 ) / Z( J4-2 )
         Z( J4 ) = Z( J4P2 )*TEMP
         DN = DNM1*TEMP
      ELSE
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) )
      end if
      DMIN = MIN( DMIN, DN )
!
      Z( J4+2 ) = DN
      Z( 4*N0-PP ) = EMIN
      RETURN
!
!     End of DLASQ6
!
      END
