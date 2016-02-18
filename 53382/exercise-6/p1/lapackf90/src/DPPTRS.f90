      SUBROUTINE DPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPPTRS solves a system of linear equations A*X = B with a symmetric
!  positive definite matrix A in packed storage using the Cholesky
!  factorization A = U**T*U or A = L*L**T computed by DPPTRF.
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
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          The triangular factor U or L from the Cholesky factorization
!          A = U**T*U or A = L*L**T, packed columnwise in a linear
!          array.  The j-th column of U or L is stored in the array AP
!          as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DTPSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( N < 0 ) then
         INFO = -2
      else if ( NRHS < 0 ) then
         INFO = -3
      else if ( LDB < MAX( 1, N ) ) then
         INFO = -6
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DPPTRS', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 .OR. NRHS.EQ.0 ) &
         RETURN
!
      if ( UPPER ) then
!
!        Solve A*X = B where A = U'*U.
!
         DO 10 I = 1, NRHS
!
!           Solve U'*X = B, overwriting B with X.
!
            CALL DTPSV( 'Upper', 'Transpose', 'Non-unit', N, AP, &
                        B( 1, I ), 1 )
!
!           Solve U*X = B, overwriting B with X.
!
            CALL DTPSV( 'Upper', 'No transpose', 'Non-unit', N, AP, &
                        B( 1, I ), 1 )
   10    CONTINUE
      ELSE
!
!        Solve A*X = B where A = L*L'.
!
         DO 20 I = 1, NRHS
!
!           Solve L*Y = B, overwriting B with X.
!
            CALL DTPSV( 'Lower', 'No transpose', 'Non-unit', N, AP, &
                        B( 1, I ), 1 )
!
!           Solve L'*X = Y, overwriting B with X.
!
            CALL DTPSV( 'Lower', 'Transpose', 'Non-unit', N, AP, &
                        B( 1, I ), 1 )
   20    CONTINUE
      end if
!
      RETURN
!
!     End of DPPTRS
!
      END
