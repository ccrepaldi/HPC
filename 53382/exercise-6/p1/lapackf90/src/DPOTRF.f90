      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
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
!  DPOTRF computes the Cholesky factorization of a real symmetric
!  positive definite matrix A.
!
!  The factorization has the form
!     A = U**T * U,  if UPLO = 'U', or
!     A = L  * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  This is the block version of the algorithm, calling Level 3 BLAS.
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
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U**T*U or A = L*L**T.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JB, NB
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DPOTF2, DSYRK, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
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
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -4
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DPOTRF', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) &
         RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 )
      if ( NB.LE.1 .OR. NB.GE.N ) then
!
!        Use unblocked code.
!
         CALL DPOTF2( UPLO, N, A, LDA, INFO )
      ELSE
!
!        Use blocked code.
!
         if ( UPPER ) then
!
!           Compute the Cholesky factorization A = U'*U.
!
            DO 10 J = 1, N, NB
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Upper', 'Transpose', JB, J-1, -ONE, &
                           A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTF2( 'Upper', JB, A( J, J ), LDA, INFO )
               if ( INFO.NE.0 ) &
                  GO TO 30
               if ( J+JB.LE.N ) then
!
!                 Compute the current block row.
!
                  CALL DGEMM( 'Transpose', 'No transpose', JB, N-J-JB+1, &
                              J-1, -ONE, A( 1, J ), LDA, A( 1, J+JB ), &
                              LDA, ONE, A( J, J+JB ), LDA )
                  CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', &
                              JB, N-J-JB+1, ONE, A( J, J ), LDA, &
                              A( J, J+JB ), LDA )
               end if
   10       CONTINUE
!
         ELSE
!
!           Compute the Cholesky factorization A = L*L'.
!
            DO 20 J = 1, N, NB
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Lower', 'No transpose', JB, J-1, -ONE, &
                           A( J, 1 ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTF2( 'Lower', JB, A( J, J ), LDA, INFO )
               if ( INFO.NE.0 ) &
                  GO TO 30
               if ( J+JB.LE.N ) then
!
!                 Compute the current block column.
!
                  CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB, &
                              J-1, -ONE, A( J+JB, 1 ), LDA, A( J, 1 ), &
                              LDA, ONE, A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'Transpose', 'Non-unit', &
                              N-J-JB+1, JB, ONE, A( J, J ), LDA, &
                              A( J+JB, J ), LDA )
               end if
   20       CONTINUE
         end if
      end if
      GO TO 40
!
   30 CONTINUE
      INFO = INFO + J - 1
!
   40 CONTINUE
      RETURN
!
!     End of DPOTRF
!
      END
