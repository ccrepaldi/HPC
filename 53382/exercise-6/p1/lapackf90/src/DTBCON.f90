      SUBROUTINE DTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, &
                         IWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            INFO, KD, LDAB, N
      DOUBLE PRECISION   RCOND
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   AB( LDAB, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DTBCON estimates the reciprocal of the condition number of a
!  triangular band matrix A, in either the 1-norm or the infinity-norm.
!
!  The norm of A is computed and an estimate is obtained for
!  norm(inv(A)), then the reciprocal of the condition number is
!  computed as
!     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies whether the 1-norm condition number or the
!          infinity-norm condition number is required:
!          = '1' or 'O':  1-norm;
!          = 'I':         Infinity-norm.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  A is upper triangular;
!          = 'L':  A is lower triangular.
!
!  DIAG    (input) CHARACTER*1
!          = 'N':  A is non-unit triangular;
!          = 'U':  A is unit triangular.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of superdiagonals or subdiagonals of the
!          triangular band matrix A.  KD >= 0.
!
!  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
!          The upper or lower triangular band matrix A, stored in the
!          first kd+1 rows of the array. The j-th column of A is stored
!          in the j-th column of the array AB as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!          If DIAG = 'U', the diagonal elements of A are not referenced
!          and are assumed to be 1.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  RCOND   (output) DOUBLE PRECISION
!          The reciprocal of the condition number of the matrix A,
!          computed as RCOND = 1/(norm(A) * norm(inv(A))).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
!
!  IWORK   (workspace) INTEGER array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOUNIT, ONENRM, UPPER
      CHARACTER          NORMIN
      INTEGER            IX, KASE, KASE1
      DOUBLE PRECISION   AINVNM, ANORM, SCALE, SMLNUM, XNORM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DLANTB
      EXTERNAL           LSAME, IDAMAX, DLAMCH, DLANTB
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLACON, DLATBS, DRSCL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      ONENRM = NORM == '1' .OR. LSAME( NORM, 'O' )
      NOUNIT = LSAME( DIAG, 'N' )
!
      if ( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) then
         INFO = -1
      else if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) then
         INFO = -2
      else if ( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) then
         INFO = -3
      else if ( N < 0 ) then
         INFO = -4
      else if ( KD < 0 ) then
         INFO = -5
      else if ( LDAB < KD+1 ) then
         INFO = -7
      end if
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DTBCON', -INFO )
         RETURN
      end if
!
!     Quick return if possible
!
      if ( N == 0 ) then
         RCOND = ONE
         RETURN
      end if
!
      RCOND = ZERO
      SMLNUM = DLAMCH( 'Safe minimum' )*DBLE( MAX( 1, N ) )
!
!     Compute the norm of the triangular matrix A.
!
      ANORM = DLANTB( NORM, UPLO, DIAG, N, KD, AB, LDAB, WORK )
!
!     Continue only if ANORM > 0.
!
      if ( ANORM.GT.ZERO ) then
!
!        Estimate the norm of the inverse of A.
!
         AINVNM = ZERO
         NORMIN = 'N'
         if ( ONENRM ) then
            KASE1 = 1
         ELSE
            KASE1 = 2
         end if
         KASE = 0
   10    CONTINUE
         CALL DLACON( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE )
         if ( KASE.NE.0 ) then
            if ( KASE == KASE1 ) then
!
!              Multiply by inv(A).
!
               CALL DLATBS( UPLO, 'No transpose', DIAG, NORMIN, N, KD, &
                            AB, LDAB, WORK, SCALE, WORK( 2*N+1 ), INFO )
            ELSE
!
!              Multiply by inv(A').
!
               CALL DLATBS( UPLO, 'Transpose', DIAG, NORMIN, N, KD, AB, &
                            LDAB, WORK, SCALE, WORK( 2*N+1 ), INFO )
            end if
            NORMIN = 'Y'
!
!           Multiply by 1/SCALE if doing so will not cause overflow.
!
            if ( SCALE.NE.ONE ) then
               IX = IDAMAX( N, WORK, 1 )
               XNORM = ABS( WORK( IX ) )
               if ( SCALE < XNORM*SMLNUM .OR. SCALE == ZERO ) &
                  GO TO 20
               CALL DRSCL( N, SCALE, WORK, 1 )
            end if
            GO TO 10
         end if
!
!        Compute the estimate of the reciprocal condition number.
!
         if ( AINVNM.NE.ZERO ) &
            RCOND = ( ONE / ANORM ) / AINVNM
      end if
!
   20 CONTINUE
      RETURN
!
!     End of DTBCON
!
      END
