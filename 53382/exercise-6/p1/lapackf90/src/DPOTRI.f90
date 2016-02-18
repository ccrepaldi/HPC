      SUBROUTINE DPOTRI( UPLO, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DPOTRI computes the inverse of a real symmetric positive definite
!  matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
!  computed by DPOTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the triangular factor U or L from the Cholesky
!          factorization A = U**T*U or A = L*L**T, as computed by
!          DPOTRF.
!          On exit, the upper or lower triangle of the (symmetric)
!          inverse of A, overwriting the input factor U or L.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the (i,i) element of the factor U or L is
!                zero, and the inverse could not be computed.
!
!  =====================================================================
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAUUM, DTRTRI, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -4
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DPOTRI', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
!
!     Invert the triangular Cholesky factor U or L.
!
      CALL DTRTRI( UPLO, 'Non-unit', N, A, LDA, INFO )
      if ( INFO.GT.0 ) &
         RETURN
!
!     Form inv(U)*inv(U)' or inv(L)'*inv(L).
!
      CALL DLAUUM( UPLO, N, A, LDA, INFO )
!
      RETURN
!
!     End of DPOTRI
!
      END
