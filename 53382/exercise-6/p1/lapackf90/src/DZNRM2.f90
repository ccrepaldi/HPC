      DOUBLE PRECISION FUNCTION DZNRM2( N, X, INCX )
!     .. Scalar Arguments ..
      INTEGER                           INCX, N
!     .. Array Arguments ..
      COMPLEX*16                        X( * )
!     ..
!
!  DZNRM2 returns the euclidean norm of a vector via the function
!  name, so that
!
!     DZNRM2 := sqrt( conjg( x' )*x )
!
!
!
!  -- This version written on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to ZLASSQ.
!     Sven Hammarling, Nag Ltd.
!
!
!     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     .. Local Scalars ..
      INTEGER               IX
      DOUBLE PRECISION      NORM, SCALE, SSQ, TEMP
!     .. Intrinsic Functions ..
      INTRINSIC             ABS, DIMAG, DBLE, SQRT
!     ..
!     .. Executable Statements ..
      if ( N < 1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE
         SCALE = ZERO
         SSQ   = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
!
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            if ( DBLE( X( IX ) ).NE.ZERO )THEN
               TEMP = ABS( DBLE( X( IX ) ) )
               if ( SCALE < TEMP )THEN
                  SSQ   = ONE   + SSQ*( SCALE/TEMP )**2
                  SCALE = TEMP
               ELSE
                  SSQ   = SSQ   +     ( TEMP/SCALE )**2
               end if
            end if
            if ( DIMAG( X( IX ) ).NE.ZERO )THEN
               TEMP = ABS( DIMAG( X( IX ) ) )
               if ( SCALE < TEMP )THEN
                  SSQ   = ONE   + SSQ*( SCALE/TEMP )**2
                  SCALE = TEMP
               ELSE
                  SSQ   = SSQ   +     ( TEMP/SCALE )**2
               end if
            end if
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      end if
!
      DZNRM2 = NORM
      RETURN
!
!     End of DZNRM2.
!
      END
