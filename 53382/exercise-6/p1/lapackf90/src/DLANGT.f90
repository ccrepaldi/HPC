      DOUBLE PRECISION FUNCTION DLANGT( NORM, N, DL, D, DU )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DL( * ), DU( * )
!     ..
!
!  Purpose
!  =======
!
!  DLANGT  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real tridiagonal matrix A.
!
!  Description
!  ===========
!
!  DLANGT returns the value
!
!     DLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANGT as described
!          above.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANGT is
!          set to zero.
!
!  DL      (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) sub-diagonal elements of A.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of A.
!
!  DU      (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) super-diagonal elements of A.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, SCALE, SUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      if ( N.LE.0 ) then
         ANORM = ZERO
      else if ( LSAME( NORM, 'M' ) ) then
!
!        Find max(abs(A(i,j))).
!
         ANORM = ABS( D( N ) )
         DO 10 I = 1, N - 1
            ANORM = MAX( ANORM, ABS( DL( I ) ) )
            ANORM = MAX( ANORM, ABS( D( I ) ) )
            ANORM = MAX( ANORM, ABS( DU( I ) ) )
   10    CONTINUE
      else if ( LSAME( NORM, 'O' ) .OR. NORM == '1' ) then
!
!        Find norm1(A).
!
         if ( N == 1 ) then
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = MAX( ABS( D( 1 ) )+ABS( DL( 1 ) ), &
                    ABS( D( N ) )+ABS( DU( N-1 ) ) )
            DO 20 I = 2, N - 1
               ANORM = MAX( ANORM, ABS( D( I ) )+ABS( DL( I ) )+ &
                       ABS( DU( I-1 ) ) )
   20       CONTINUE
         end if
      else if ( LSAME( NORM, 'I' ) ) then
!
!        Find normI(A).
!
         if ( N == 1 ) then
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = MAX( ABS( D( 1 ) )+ABS( DU( 1 ) ), &
                    ABS( D( N ) )+ABS( DL( N-1 ) ) )
            DO 30 I = 2, N - 1
               ANORM = MAX( ANORM, ABS( D( I ) )+ABS( DU( I ) )+ &
                       ABS( DL( I-1 ) ) )
   30       CONTINUE
         end if
      else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) then
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         CALL DLASSQ( N, D, 1, SCALE, SUM )
         if ( N.GT.1 ) then
            CALL DLASSQ( N-1, DL, 1, SCALE, SUM )
            CALL DLASSQ( N-1, DU, 1, SCALE, SUM )
         end if
         ANORM = SCALE*SQRT( SUM )
      end if
!
      DLANGT = ANORM
      RETURN
!
!     End of DLANGT
!
      END
