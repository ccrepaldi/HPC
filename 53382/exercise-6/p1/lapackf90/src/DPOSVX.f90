      SUBROUTINE DPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, &
                         S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, &
                         IWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          EQUED, FACT, UPLO
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
      DOUBLE PRECISION   RCOND
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ), &
                         BERR( * ), FERR( * ), S( * ), WORK( * ), &
                         X( LDX, * )
!     ..
!
!  Purpose
!  =======
!
!  DPOSVX uses the Cholesky factorization A = U**T*U or A = L*L**T to
!  compute the solution to a real system of linear equations
!     A * X = B,
!  where A is an N-by-N symmetric positive definite matrix and X and B
!  are N-by-NRHS matrices.
!
!  Error bounds on the solution and a condition estimate are also
!  provided.
!
!  Description
!  ===========
!
!  The following steps are performed:
!
!  1. If FACT = 'E', real scaling factors are computed to equilibrate
!     the system:
!        diag(S) * A * diag(S) * inv(diag(S)) * X = diag(S) * B
!     Whether or not the system will be equilibrated depends on the
!     scaling of the matrix A, but if equilibration is used, A is
!     overwritten by diag(S)*A*diag(S) and B by diag(S)*B.
!
!  2. If FACT = 'N' or 'E', the Cholesky decomposition is used to
!     factor the matrix A (after equilibration if FACT = 'E') as
!        A = U**T* U,  if UPLO = 'U', or
!        A = L * L**T,  if UPLO = 'L',
!     where U is an upper triangular matrix and L is a lower triangular
!     matrix.
!
!  3. If the leading i-by-i principal minor is not positive definite,
!     then the routine returns with INFO = i. Otherwise, the factored
!     form of A is used to estimate the condition number of the matrix
!     A.  If the reciprocal of the condition number is less than machine
!     precision, INFO = N+1 is returned as a warning, but the routine
!     still goes on to solve for X and compute error bounds as
!     described below.
!
!  4. The system of equations is solved for X using the factored form
!     of A.
!
!  5. Iterative refinement is applied to improve the computed solution
!     matrix and calculate error bounds and backward error estimates
!     for it.
!
!  6. If equilibration was used, the matrix X is premultiplied by
!     diag(S) so that it solves the original system before
!     equilibration.
!
!  Arguments
!  =========
!
!  FACT    (input) CHARACTER*1
!          Specifies whether or not the factored form of the matrix A is
!          supplied on entry, and if not, whether the matrix A should be
!          equilibrated before it is factored.
!          = 'F':  On entry, AF contains the factored form of A.
!                  If EQUED = 'Y', the matrix A has been equilibrated
!                  with scaling factors given by S.  A and AF will not
!                  be modified.
!          = 'N':  The matrix A will be copied to AF and factored.
!          = 'E':  The matrix A will be equilibrated if necessary, then
!                  copied to AF and factored.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrices B and X.  NRHS >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A, except if FACT = 'F' and
!          EQUED = 'Y', then A must contain the equilibrated matrix
!          diag(S)*A*diag(S).  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.  A is not modified if
!          FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit.
!
!          On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by
!          diag(S)*A*diag(S).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  AF      (input or output) DOUBLE PRECISION array, dimension (LDAF,N)
!          If FACT = 'F', then AF is an input argument and on entry
!          contains the triangular factor U or L from the Cholesky
!          factorization A = U**T*U or A = L*L**T, in the same storage
!          format as A.  If EQUED .ne. 'N', then AF is the factored form
!          of the equilibrated matrix diag(S)*A*diag(S).
!
!          If FACT = 'N', then AF is an output argument and on exit
!          returns the triangular factor U or L from the Cholesky
!          factorization A = U**T*U or A = L*L**T of the original
!          matrix A.
!
!          If FACT = 'E', then AF is an output argument and on exit
!          returns the triangular factor U or L from the Cholesky
!          factorization A = U**T*U or A = L*L**T of the equilibrated
!          matrix A (see the description of A for the form of the
!          equilibrated matrix).
!
!  LDAF    (input) INTEGER
!          The leading dimension of the array AF.  LDAF >= max(1,N).
!
!  EQUED   (input or output) CHARACTER*1
!          Specifies the form of equilibration that was done.
!          = 'N':  No equilibration (always true if FACT = 'N').
!          = 'Y':  Equilibration was done, i.e., A has been replaced by
!                  diag(S) * A * diag(S).
!          EQUED is an input argument if FACT = 'F'; otherwise, it is an
!          output argument.
!
!  S       (input or output) DOUBLE PRECISION array, dimension (N)
!          The scale factors for A; not accessed if EQUED = 'N'.  S is
!          an input argument if FACT = 'F'; otherwise, S is an output
!          argument.  If FACT = 'F' and EQUED = 'Y', each element of S
!          must be positive.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS right hand side matrix B.
!          On exit, if EQUED = 'N', B is not modified; if EQUED = 'Y',
!          B is overwritten by diag(S) * B.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)
!          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X to
!          the original system of equations.  Note that if EQUED = 'Y',
!          A and B are modified on exit, and the solution to the
!          equilibrated system is inv(diag(S))*X.
!
!  LDX     (input) INTEGER
!          The leading dimension of the array X.  LDX >= max(1,N).
!
!  RCOND   (output) DOUBLE PRECISION
!          The estimate of the reciprocal condition number of the matrix
!          A after equilibration (if done).  If RCOND is less than the
!          machine precision (in particular, if RCOND = 0), the matrix
!          is singular to working precision.  This condition is
!          indicated by a return code of INFO > 0.
!
!  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)
!          The estimated forward error bound for each solution vector
!          X(j) (the j-th column of the solution matrix X).
!          If XTRUE is the true solution corresponding to X(j), FERR(j)
!          is an estimated upper bound for the magnitude of the largest
!          element in (X(j) - XTRUE) divided by the magnitude of the
!          largest element in X(j).  The estimate is as reliable as
!          the estimate for RCOND, and is almost always a slight
!          overestimate of the true error.
!
!  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)
!          The componentwise relative backward error of each solution
!          vector X(j) (i.e., the smallest relative change in
!          any element of A or B that makes X(j) an exact solution).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
!
!  IWORK   (workspace) INTEGER array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, and i is
!                <= N:  the leading minor of order i of A is
!                       not positive definite, so the factorization
!                       could not be completed, and the solution has not
!                       been computed. RCOND = 0 is returned.
!                = N+1: U is nonsingular, but RCOND is less than machine
!                       precision, meaning that the matrix is singular
!                       to working precision.  Nevertheless, the
!                       solution and error bounds are computed because
!                       there are a number of situations where the
!                       computed solution can be more accurate than the
!                       value of RCOND would suggest.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            EQUIL, NOFACT, RCEQU
      INTEGER            I, INFEQU, J
      DOUBLE PRECISION   AMAX, ANORM, BIGNUM, SCOND, SMAX, SMIN, SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           LSAME, DLAMCH, DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLACPY, DLAQSY, DPOCON, DPOEQU, DPORFS, DPOTRF, &
                         DPOTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      if ( NOFACT .OR. EQUIL ) then
         EQUED = 'N'
         RCEQU = .FALSE.
      ELSE
         RCEQU = LSAME( EQUED, 'Y' )
         SMLNUM = DLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
      end if
!
!     Test the input parameters.
!
      if ( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT.LSAME( FACT, 'F' ) ) &
           THEN
         INFO = -1
      else if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) &
                THEN
         INFO = -2
      else if ( N < 0 ) then
         INFO = -3
      else if ( NRHS < 0 ) then
         INFO = -4
      else if ( LDA < MAX( 1, N ) ) then
         INFO = -6
      else if ( LDAF < MAX( 1, N ) ) then
         INFO = -8
      else if ( LSAME( FACT, 'F' ) .AND. .NOT. &
               ( RCEQU .OR. LSAME( EQUED, 'N' ) ) ) then
         INFO = -9
      ELSE
         if ( RCEQU ) then
            SMIN = BIGNUM
            SMAX = ZERO
            DO 10 J = 1, N
               SMIN = MIN( SMIN, S( J ) )
               SMAX = MAX( SMAX, S( J ) )
   10       CONTINUE
            if ( SMIN.LE.ZERO ) then
               INFO = -10
            else if ( N.GT.0 ) then
               SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )
            ELSE
               SCOND = ONE
            end if
         end if
         if ( INFO == 0 ) then
            if ( LDB < MAX( 1, N ) ) then
               INFO = -12
            else if ( LDX < MAX( 1, N ) ) then
               INFO = -14
            end if
         end if
      end if
!
      if ( INFO.NE.0 ) then
         CALL XERBLA( 'DPOSVX', -INFO )
         RETURN
      end if
!
      if ( EQUIL ) then
!
!        Compute row and column scalings to equilibrate the matrix A.
!
         CALL DPOEQU( N, A, LDA, S, SCOND, AMAX, INFEQU )
         if ( INFEQU == 0 ) then
!
!           Equilibrate the matrix.
!
            CALL DLAQSY( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )
            RCEQU = LSAME( EQUED, 'Y' )
         end if
      end if
!
!     Scale the right hand side.
!
      if ( RCEQU ) then
         DO 30 J = 1, NRHS
            DO 20 I = 1, N
               B( I, J ) = S( I )*B( I, J )
   20       CONTINUE
   30    CONTINUE
      end if
!
      if ( NOFACT .OR. EQUIL ) then
!
!        Compute the Cholesky factorization A = U'*U or A = L*L'.
!
         CALL DLACPY( UPLO, N, N, A, LDA, AF, LDAF )
         CALL DPOTRF( UPLO, N, AF, LDAF, INFO )
!
!        Return if INFO is non-zero.
!
         if ( INFO.NE.0 ) then
            if ( INFO.GT.0 ) &
               RCOND = ZERO
            RETURN
         end if
      end if
!
!     Compute the norm of the matrix A.
!
      ANORM = DLANSY( '1', UPLO, N, A, LDA, WORK )
!
!     Compute the reciprocal of the condition number of A.
!
      CALL DPOCON( UPLO, N, AF, LDAF, ANORM, RCOND, WORK, IWORK, INFO )
!
!     Set INFO = N+1 if the matrix is singular to working precision.
!
      if ( RCOND < DLAMCH( 'Epsilon' ) ) &
         INFO = N + 1
!
!     Compute the solution matrix X.
!
      CALL DLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL DPOTRS( UPLO, N, NRHS, AF, LDAF, X, LDX, INFO )
!
!     Use iterative refinement to improve the computed solution and
!     compute error bounds and backward error estimates for it.
!
      CALL DPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, LDX, &
                   FERR, BERR, WORK, IWORK, INFO )
!
!     Transform the solution matrix X to a solution of the original
!     system.
!
      if ( RCEQU ) then
         DO 50 J = 1, NRHS
            DO 40 I = 1, N
               X( I, J ) = S( I )*X( I, J )
   40       CONTINUE
   50    CONTINUE
         DO 60 J = 1, NRHS
            FERR( J ) = FERR( J ) / SCOND
   60    CONTINUE
      end if
!
      RETURN
!
!     End of DPOSVX
!
      END
