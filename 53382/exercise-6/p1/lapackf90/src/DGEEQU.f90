      SUBROUTINE DGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, &
                         INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
      DOUBLE PRECISION   AMAX, COLCND, ROWCND
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), R( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEEQU computes row and column scalings intended to equilibrate an
!  M-by-N matrix A and reduce its condition number.  R returns the row
!  scale factors and C the column scale factors, chosen to try to make
!  the largest element in each row and column of the matrix B with
!  elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.
!
!  R(i) and C(j) are restricted to be between SMLNUM = smallest safe
!  number and BIGNUM = largest safe number.  Use of these scaling
!  factors is not guaranteed to reduce the condition number of A but
!  works well in practice.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The M-by-N matrix whose equilibration factors are
!          to be computed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  R       (output) DOUBLE PRECISION array, dimension (M)
!          If INFO = 0 or INFO > M, R contains the row scale factors
!          for A.
!
!  C       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0,  C contains the column scale factors for A.
!
!  ROWCND  (output) DOUBLE PRECISION
!          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
!          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
!          AMAX is neither too large nor too small, it is not worth
!          scaling by R.
!
!  COLCND  (output) DOUBLE PRECISION
!          If INFO = 0, COLCND contains the ratio of the smallest
!          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
!          worth scaling by C.
!
!  AMAX    (output) DOUBLE PRECISION
!          Absolute value of largest matrix element.  If AMAX is very
!          close to overflow or very close to underflow, the matrix
!          should be scaled.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i,  and i is
!                <= M:  the i-th row of A is exactly zero
!                >  M:  the (i-M)-th column of A is exactly zero
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   BIGNUM, RCMAX, RCMIN, SMLNUM
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      if ( M < 0 ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDA < MAX( 1, M ) ) then
         INFO = -4
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DGEEQU', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( M == 0 .OR. N.EQ.0 ) then
         ROWCND = ONE
         COLCND = ONE
         AMAX = ZERO
         RETURN
      end if
!
!     Get machine constants.
!
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
!
!     Compute row scale factors.
!
      R(1:m) = ZERO
!
!     Find the maximum element in each row.
!
      DO 30 J = 1, N
         DO 20 I = 1, M
            R( I ) = MAX( R( I ), ABS( A( I, J ) ) )
   20    CONTINUE
   30 CONTINUE
!
!     Find the maximum and minimum scale factors.
!
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 40 I = 1, M
         RCMAX = MAX( RCMAX, R( I ) )
         RCMIN = MIN( RCMIN, R( I ) )
   40 CONTINUE
      AMAX = RCMAX
!
      if ( RCMIN == ZERO ) then
!
!        Find the first zero scale factor and return an error code.
!
         DO 50 I = 1, M
            if ( R( I ) == ZERO ) then
               INFO = I
               RETURN
            end if
   50    CONTINUE
      ELSE
!
!        Invert the scale factors.
!
         DO 60 I = 1, M
            R( I ) = ONE / MIN( MAX( R( I ), SMLNUM ), BIGNUM )
   60    CONTINUE
!
!        Compute ROWCND = min(R(I)) / max(R(I))
!
         ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      end if
!
!     Compute column scale factors
!
      DO 70 J = 1, N
         C( J ) = ZERO
   70 CONTINUE
!
!     Find the maximum element in each column,
!     assuming the row scaling computed above.
!
      DO 90 J = 1, N
         DO 80 I = 1, M
            C( J ) = MAX( C( J ), ABS( A( I, J ) )*R( I ) )
   80    CONTINUE
   90 CONTINUE
!
!     Find the maximum and minimum scale factors.
!
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO J = 1, N
         RCMIN = MIN( RCMIN, C( J ) )
         RCMAX = MAX( RCMAX, C( J ) )
      end do

      if ( RCMIN == ZERO ) then
!
!        Find the first zero scale factor and return an error code.
!
         DO 110 J = 1, N
            if ( C( J ) == ZERO ) then
               INFO = M + J
               RETURN
            end if
  110    CONTINUE
      ELSE
!
!        Invert the scale factors.
!
         DO 120 J = 1, N
            C( J ) = ONE / MIN( MAX( C( J ), SMLNUM ), BIGNUM )
  120    CONTINUE
!
!        Compute COLCND = min(C(J)) / max(C(J))
!
         COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      end if
!
      RETURN
!
!     End of DGEEQU
!
      END
