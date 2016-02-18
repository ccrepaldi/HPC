      SUBROUTINE DPOEQU( N, A, LDA, S, SCOND, AMAX, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   AMAX, SCOND
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), S( * )
!     ..
!
!  Purpose
!  =======
!
!  DPOEQU computes row and column scalings intended to equilibrate a
!  symmetric positive definite matrix A and reduce its condition number
!  (with respect to the two-norm).  S contains the scale factors,
!  S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with
!  elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This
!  choice of S puts the condition number of B within a factor N of the
!  smallest possible condition number over all possible diagonal
!  scalings.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The N-by-N symmetric positive definite matrix whose scaling
!          factors are to be computed.  Only the diagonal elements of A
!          are referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  S       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, S contains the scale factors for A.
!
!  SCOND   (output) DOUBLE PRECISION
!          If INFO = 0, S contains the ratio of the smallest S(i) to
!          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too
!          large nor too small, it is not worth scaling by S.
!
!  AMAX    (output) DOUBLE PRECISION
!          Absolute value of largest matrix element.  If AMAX is very
!          close to overflow or very close to underflow, the matrix
!          should be scaled.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the i-th diagonal element is nonpositive.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   SMIN
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      if ( N < 0 ) then
         INFO = -1
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -3
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DPOEQU', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) then
         SCOND = ONE
         AMAX = ZERO
         RETURN
      end if
!
!     Find the minimum and maximum diagonal elements.
!
      S( 1 ) = A( 1, 1 )
      SMIN = S( 1 )
      AMAX = S( 1 )
      DO 10 I = 2, N
         S( I ) = A( I, I )
         SMIN = MIN( SMIN, S( I ) )
         AMAX = MAX( AMAX, S( I ) )
   10 CONTINUE
!
      if ( SMIN.LE.ZERO ) then
!
!        Find the first non-positive diagonal element and return.
!
         DO 20 I = 1, N
            if ( S( I ).LE.ZERO ) then
               INFO = I
               RETURN
            end if
   20    CONTINUE
      ELSE
!
!        Set the scale factors to the reciprocals
!        of the diagonal elements.
!
         DO 30 I = 1, N
            S( I ) = ONE / SQRT( S( I ) )
   30    CONTINUE
!
!        Compute SCOND = min(S(I)) / max(S(I))
!
         SCOND = SQRT( SMIN ) / SQRT( AMAX )
      end if
      RETURN
!
!     End of DPOEQU
!
      END
