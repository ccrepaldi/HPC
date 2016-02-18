      DOUBLE PRECISION FUNCTION DLANSB( NORM, UPLO, N, K, AB, LDAB, &
                       WORK )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            K, LDAB, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLANSB  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the element of  largest absolute value  of an
!  n by n symmetric band matrix A,  with k super-diagonals.
!
!  Description
!  ===========
!
!  DLANSB returns the value
!
!     DLANSB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!          Specifies the value to be returned in DLANSB as described
!          above.
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          band matrix A is supplied.
!          = 'U':  Upper triangular part is supplied
!          = 'L':  Lower triangular part is supplied
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANSB is
!          set to zero.
!
!  K       (input) INTEGER
!          The number of super-diagonals or sub-diagonals of the
!          band matrix A.  K >= 0.
!
!  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
!          The upper or lower triangle of the symmetric band matrix A,
!          stored in the first K+1 rows of AB.  The j-th column of A is
!          stored in the j-th column of the array AB as follows:
!          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= K+1.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!          WORK is not referenced.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, L
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASSQ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      if ( N == 0 ) then
         VALUE = ZERO
      else if ( LSAME( NORM, 'M' ) ) then
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         if ( LSAME( UPLO, 'U' ) ) then
            DO 20 J = 1, N
               DO 10 I = MAX( K+2-J, 1 ), K + 1
                  VALUE = MAX( VALUE, ABS( AB( I, J ) ) )
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               DO 30 I = 1, MIN( N+1-J, K+1 )
                  VALUE = MAX( VALUE, ABS( AB( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         end if
      else if ( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. &
               ( NORM == '1' ) ) then
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
         VALUE = ZERO
         if ( LSAME( UPLO, 'U' ) ) then
            DO 60 J = 1, N
               SUM = ZERO
               L = K + 1 - J
               DO 50 I = MAX( 1, J-K ), J - 1
                  ABSA = ABS( AB( L+I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( AB( K+1, J ) )
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( AB( 1, J ) )
               L = 1 - J
               DO 90 I = J + 1, MIN( N, J+K )
                  ABSA = ABS( AB( L+I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         end if
      else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) then
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         if ( K.GT.0 ) then
            if ( LSAME( UPLO, 'U' ) ) then
               DO 110 J = 2, N
                  CALL DLASSQ( MIN( J-1, K ), AB( MAX( K+2-J, 1 ), J ), &
                               1, SCALE, SUM )
  110          CONTINUE
               L = K + 1
            ELSE
               DO 120 J = 1, N - 1
                  CALL DLASSQ( MIN( N-J, K ), AB( 2, J ), 1, SCALE, &
                               SUM )
  120          CONTINUE
               L = 1
            end if
            SUM = 2*SUM
         ELSE
            L = 1
         end if
         CALL DLASSQ( N, AB( L, 1 ), LDAB, SCALE, SUM )
         VALUE = SCALE*SQRT( SUM )
      end if
!
      DLANSB = VALUE
      RETURN
!
!     End of DLANSB
!
      END
