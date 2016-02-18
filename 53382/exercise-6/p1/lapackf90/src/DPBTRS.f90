      SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPBTRS solves a system of linear equations A*X = B with a symmetric
!  positive definite band matrix A using the Cholesky factorization
!  A = U**T*U or A = L*L**T computed by DPBTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangular factor stored in AB;
!          = 'L':  Lower triangular factor stored in AB.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
!          The triangular factor U or L from the Cholesky factorization
!          A = U**T*U or A = L*L**T of the band matrix A, stored in the
!          first KD+1 rows of the array.  The j-th column of U or L is
!          stored in the j-th column of the array AB as follows:
!          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
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
      INTEGER            J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DTBSV, XERBLA
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
      else if ( KD < 0 ) then
         INFO = -3
      else if ( NRHS < 0 ) then
         INFO = -4
      else if ( LDAB < KD+1 ) then
         INFO = -6
      else if ( LDB < MAX( 1, N ) ) then
         INFO = -8
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DPBTRS', -INFO )
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
         DO 10 J = 1, NRHS
!
!           Solve U'*X = B, overwriting B with X.
!
            CALL DTBSV( 'Upper', 'Transpose', 'Non-unit', N, KD, AB, &
                        LDAB, B( 1, J ), 1 )
!
!           Solve U*X = B, overwriting B with X.
!
            CALL DTBSV( 'Upper', 'No transpose', 'Non-unit', N, KD, AB, &
                        LDAB, B( 1, J ), 1 )
   10    CONTINUE
      ELSE
!
!        Solve A*X = B where A = L*L'.
!
         DO 20 J = 1, NRHS
!
!           Solve L*X = B, overwriting B with X.
!
            CALL DTBSV( 'Lower', 'No transpose', 'Non-unit', N, KD, AB, &
                        LDAB, B( 1, J ), 1 )
!
!           Solve L'*X = B, overwriting B with X.
!
            CALL DTBSV( 'Lower', 'Transpose', 'Non-unit', N, KD, AB, &
                        LDAB, B( 1, J ), 1 )
   20    CONTINUE
      end if
!
      RETURN
!
!     End of DPBTRS
!
      END
