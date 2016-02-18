      SUBROUTINE DTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DTPTRS solves a triangular system of the form
!
!     A * X = B  or  A**T * X = B,
!
!  where A is a triangular matrix of order N stored in packed format,
!  and B is an N-by-NRHS matrix.  A check is made to verify that A is
!  nonsingular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  A is upper triangular;
!          = 'L':  A is lower triangular.
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A**T * X = B  (Transpose)
!          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
!
!  DIAG    (input) CHARACTER*1
!          = 'N':  A is non-unit triangular;
!          = 'U':  A is unit triangular.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          The upper or lower triangular matrix A, packed columnwise in
!          a linear array.  The j-th column of A is stored in the array
!          AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, if INFO = 0, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the i-th diagonal element of A is zero,
!                indicating that the matrix is singular and the
!                solutions X have not been computed.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J, JC
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
      NOUNIT = LSAME( DIAG, 'N' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -1
      else if ( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT. &
               LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) then
         INFO = -2
      else if ( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) then
         INFO = -3
      else if ( N < 0 ) then
         INFO = -4
      else if ( NRHS < 0 ) then
         INFO = -5
      else if ( LDB < MAX( 1, N ) ) then
         INFO = -8
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DTPTRS', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
!
!     Check for singularity.
!
      if ( NOUNIT ) then
         if ( UPPER ) then
            JC = 1
            DO 10 INFO = 1, N
               if ( AP( JC+INFO-1 ) == ZERO ) &
                  RETURN
               JC = JC + INFO
   10       CONTINUE
         ELSE
            JC = 1
            DO 20 INFO = 1, N
               if ( AP( JC ) == ZERO ) &
                  RETURN
               JC = JC + N - INFO + 1
   20       CONTINUE
         end if
      end if
      INFO = 0
!
!     Solve A * x = b  or  A' * x = b.
!
      DO 30 J = 1, NRHS
         CALL DTPSV( UPLO, TRANS, DIAG, N, AP, B( 1, J ), 1 )
   30 CONTINUE
!
      RETURN
!
!     End of DTPTRS
!
      END
