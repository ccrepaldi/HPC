      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1998
!
!     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  INTEGER
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, &
                         NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
      IEEECK = 1
!
      POSINF = ONE / ZERO
      if ( POSINF.LE.ONE ) then
         IEEECK = 0
         RETURN
      end if
!
      NEGINF = -ONE / ZERO
      if ( NEGINF.GE.ZERO ) then
         IEEECK = 0
         RETURN
      end if
!
      NEGZRO = ONE / ( NEGINF+ONE )
      if ( NEGZRO.NE.ZERO ) then
         IEEECK = 0
         RETURN
      end if
!
      NEGINF = ONE / NEGZRO
      if ( NEGINF.GE.ZERO ) then
         IEEECK = 0
         RETURN
      end if
!
      NEWZRO = NEGZRO + ZERO
      if ( NEWZRO.NE.ZERO ) then
         IEEECK = 0
         RETURN
      end if
!
      POSINF = ONE / NEWZRO
      if ( POSINF.LE.ONE ) then
         IEEECK = 0
         RETURN
      end if
!
      NEGINF = NEGINF*POSINF
      if ( NEGINF.GE.ZERO ) then
         IEEECK = 0
         RETURN
      end if
!
      POSINF = POSINF*POSINF
      if ( POSINF.LE.ONE ) then
         IEEECK = 0
         RETURN
      end if
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      if ( ISPEC == 0 ) &
         RETURN
!
      NAN1 = POSINF + NEGINF
!
      NAN2 = POSINF / NEGINF
!
      NAN3 = POSINF / POSINF
!
      NAN4 = POSINF*ZERO
!
      NAN5 = NEGINF*NEGZRO
!
      NAN6 = NAN5*0.0
!
      if ( NAN1 == NAN1 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN2 == NAN2 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN3 == NAN3 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN4 == NAN4 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN5 == NAN5 ) then
         IEEECK = 0
         RETURN
      end if
!
      if ( NAN6 == NAN6 ) then
         IEEECK = 0
         RETURN
      end if
!
      RETURN
      END
